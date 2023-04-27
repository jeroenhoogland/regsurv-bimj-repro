# This is the main simulation script

# Within master_sim.R, the following should be specified in the Global environment 
# `iter` (the iteration number) 
# `sample.size.settings` (vector with elements from 1:4 for the sample size settings to be 
#   evaluated (N=100, 250, 500, 1000 respectively)
# `highercens` (FALSE for main simulations and true for increased censoring rate)
# `higherp` (FALSE for main simulations and true for increased number of covariates)

library(dplyr)
library(regsurv)
library(CVXR)
library(simsurv)
library(survival)
library(MASS)
library(glmnet)
library(rstpm2)
library(polspline)
library(Hmisc)

load(file="data/groundtruth.RData")
load("data/validation.RData")
source("sim/helperFunctions.R")

n <- c(100, 250, 500, 1000)
nfolds <- 10
maxit <- 500L
nvarobs <- 10
rmse <- e90 <- array(NA, dim=c(4, 10, 7))
time.per.model <- array(NA, dim=c(4, 10))
cal <- cstats <- brier <- cstatstd <- cstats.fixed <- lapply(1:length(n), function(x){vector(mode = "list", length = 10)})
coxzph.results <- list()
time.elapsed.sec <- rep(NA, length(n))
fixed.times <- c(2.5, 5, 7.5, 10, 20, 30)
dimnames(rmse)  <- dimnames(e90) <- list(
  paste0("n", n),
  c("logh.logt.pen", "logH.logt.pen", "coxph", "coxph.tt", "cox.lasso",
    "cox.ridge.tt", "stpm2.ph", "stpm2.nonph", "pstpm2.ph", "pstpm2.nonph"),
  c("overall", paste0("t", fixed.times)))

## Pre-compute the various data representations that are required for all prediction models
if(higherp){
  pred.data.regsurv <- list(as.matrix(data.frame("eventtime"=validation$eventtime,
                                                 validation[ ,c(grep("x", names(validation))[1:nvarobs],
                                                                grep("n[0-9]", names(validation)))])))
  pred.data.regsurv[2:7] <- lapply(1:length(fixed.times), function(t){
    as.matrix(data.frame("eventtime"=rep(fixed.times[t], length(validation$eventtime)),
                         validation[ ,c(grep("x", names(validation))[1:nvarobs],
                                        grep("n[0-9]", names(validation)))]))})
  
  pred.data.other <- list(validation)
  pred.data.other[2:7] <- lapply(1:length(fixed.times), function(t){
    val.fixed <- validation
    val.fixed$eventtime <- fixed.times[t]
    val.fixed})
} else {
  pred.data.regsurv <- list(as.matrix(data.frame("eventtime"=validation$eventtime,
                                                 validation[ ,grep("x", names(validation))[1:nvarobs]])))
  pred.data.regsurv[2:7] <- lapply(1:length(fixed.times), function(t){
    as.matrix(data.frame("eventtime"=rep(fixed.times[t], length(validation$eventtime)),
                         validation[ ,grep("x", names(validation))[1:nvarobs]]))})
  
  pred.data.other <- list(validation[ -grep("n[0-9]", names(validation))])
  pred.data.other[2:7] <- lapply(1:length(fixed.times), function(t){
    val.fixed <- validation[ -grep("n[0-9]", names(validation))]
    val.fixed$eventtime <- fixed.times[t]
    val.fixed})
}

## Store the ground truth survival probabilities in Strue for easy access
Strue <- cbind(exp(-groundtruth$Htrue), exp(-groundtruth$Htrue2.5), exp(-groundtruth$Htrue5), exp(-groundtruth$Htrue7.5),
               exp(-groundtruth$Htrue10), exp(-groundtruth$Htrue20), exp(-groundtruth$Htrue30))

## a grid of all evaluation times
time.grid <- cbind(validation$eventtime, fixed.times[1], fixed.times[2], fixed.times[3], fixed.times[4], fixed.times[5], fixed.times[6])

# grid for the evaluation of calibration performance (for each time point)
predict.grid <- sapply(1:length(fixed.times), function(t) quantile(1-Strue[ ,t+1],probs=seq(0.01,.99,.01)))

s=1

# object names for those objects that should remain in memory
ws.names <- c(ls(), "time.elapsed.sec", "ws.names")

for(s in sample.size.settings){

  ptm <- proc.time()
  
  # sample development data of size n[s] from the dev_population
  load("data/dev_population.RData")
  simdata <- dev_population[sample(1:nrow(dev_population), n[s]), ]
  
  if(!higherp){
    simdata <- simdata[ -grep("n[0-9]", names(simdata))]
  }
  
  # add exponential censoring to get 50% censoring in the data set
  if(highercens){
    C <- rexp(n=nrow(simdata), rate=.067) # 0.067 chosen to get ~50%
    simdata$status <- as.numeric(simdata$eventtime < C & simdata$status==1)
    simdata$eventtime <- pmin(simdata$eventtime, C)
  }
  
  rm(dev_population) # avoid having the entire data set in memory all the time

  ntimebasis <- 5
  nitimebasis <- 2

  ###################################
  ## log time log hazard           ##
  ###################################

  model.time <- proc.time()
  
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,c(grep("x", names(simdata))[1:nvarobs],
                                         grep("n[0-9]", names(simdata)))]),
                   model.scale="loghazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=ntimebasis,
                   tv=1:nvarobs,
                   nitimebasis=nitimebasis,
                   qpoints=25)

  groups <- list()
  groups[[1]] <- 1
  groups[[2]] <- 2:(ntimebasis+1)
  for(i in 0:(length(prep$which.param[[2]])-1)){
    groups[[3+i]] <- prep$which.param[[2]][i+1]
  }
  a <- 1:2
  for(i in 0:(nvarobs-1)){
    groups[[3+(length(prep$which.param[[2]]))+i]] <- prep$which.param[[3]][a+(2*i)]
  }

  penpars <- c(0,rep(1, length(groups)-1))
  l1l2 <- c(rep(0, length(groups)-10-as.numeric(higherp)*20),rep(1, 10+as.numeric(higherp)*20))

  obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds,
                       maxit=maxit, pred.data=pred.data, lambda.grid = exp(seq(-2,4,.5)))

  if(is.null(obj)){
    obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds,
                         maxit=maxit, pred.data=pred.data, lambda.grid = exp(seq(-2,4,.5)))
  }

  time.per.model[s,1] <- (proc.time() - model.time)["elapsed"]

  if(!is.null(obj)){
    Shat <- sapply(1:7, function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                            newdata=pred.data.regsurv[[x]], type="surv"))
    rmse[s, 1, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
    e90[s, 1, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
    cal[[s]][[1]] <-  tryCatch(
      expr=lapply(1:length(fixed.times), function(t)
        calibrate(Shat=Shat[ ,1+t],
                  time=fixed.times[t], validation=validation,
                  predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
      error=function(e) NULL)
    cstats[[s]][[1]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
    brier[[s]][[1]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
    cstatstd[[s]][[1]] <- mean(sapply(1:10, function(x){
      ss <- sample(1:nrow(validation), 100)
      ctd(validation$eventtime[ss], validation$status[ss],
          X=pred.data.regsurv[[1]][ss ,-which(colnames(pred.data.regsurv[[1]]) == "eventtime")],
          obj$mod, type="regsurv", prep=prep, lambda.index=obj$cv$lambda.min.index)$Ctd
    }))
    cstats.fixed[[s]][[1]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

    rm(Shat)
    rm(obj)
  }

  ###################################
  ## log time log Hazard           ##
  ###################################

  model.time <- proc.time()

  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,c(grep("x", names(simdata))[1:nvarobs],
                                           grep("n[0-9]", names(simdata)))]),
                   model.scale="logHazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=ntimebasis,
                   tv=1:nvarobs,
                   nitimebasis=nitimebasis,
                   qpoints=NULL)

  obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds,
                       maxit=maxit, pred.data=pred.data, lambda.grid = exp(seq(-2,4,.5)))

  if(is.null(obj)){
    obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds,
                         maxit=maxit, pred.data=pred.data, lambda.grid = exp(seq(-2,4,.5)))
  }

  time.per.model[s,2] <- (proc.time() - model.time)["elapsed"]

  if(!is.null(obj)){
    Shat <- sapply(1:7, function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                            newdata=pred.data.regsurv[[x]], type="surv"))
    rmse[s, 2, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
    e90[s, 2, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
    cal[[s]][[2]] <- tryCatch(
      expr=lapply(1:length(fixed.times), function(t)
        calibrate(Shat=Shat[ ,1+t],
                  time=fixed.times[t], validation=validation,
                  predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
      error=function(e) NULL)
    cstats[[s]][[2]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
    brier[[s]][[2]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
    cstatstd[[s]][[2]] <- mean(sapply(1:10, function(x){
      ss <- sample(1:nrow(validation), 100)
      ctd(validation$eventtime[ss], validation$status[ss],
          X=pred.data.regsurv[[1]][ss ,-which(colnames(pred.data.regsurv[[1]]) == "eventtime")],
          obj$mod, type="regsurv", prep=prep, lambda.index=obj$cv$lambda.min.index)$Ctd
    }))
    cstats.fixed[[s]][[2]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

    rm(Shat)
    rm(obj)
  }

  ###################################
  ## coxph                         ##
  ###################################

  model.time <- proc.time()

  if(!higherp){
    form <- formula(paste("Surv(eventtime, status) ~ ", paste0("x", 1:nvarobs, collapse="+")))
  } else {
    form <- formula(paste("Surv(eventtime, status) ~ ", paste0("x", 1:nvarobs, collapse="+"),
                          "+", paste0("n", 1:20, collapse="+")))
  }

  # see coxph.control() help for timefix (especially in simulated data without ties)
  mod.coxph <- coxph(form, data=simdata, control=coxph.control(timefix=FALSE))

  time.per.model[s,3] <- (proc.time() - model.time)["elapsed"]

  coxzph.results[[s]] <- cox.zph(mod.coxph)$table

  Shat <- sapply(1:7, function(x){
    exp(-predict(mod.coxph, newdata=pred.data.other[[x]], type="expected"))
  })

  rmse[s, 3, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 3, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[3]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[3]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[3]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[3]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        mod.coxph, type="coxph")$Ctd
  }))
  cstats.fixed[[s]][[3]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat)

  ###################################
  ## coxph.tt                      ##
  ###################################

  model.time <- proc.time()

  # https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
  dtimes <- unique(quantile(simdata$eventtime[simdata$status==1], probs=seq(0,1,length.out=101), type=1))
  simdata.tt <- survSplit(Surv(eventtime, status==1) ~., simdata, cut=dtimes)
  names(simdata.tt)[names(simdata.tt) == "event"] <- "status"
  time.spline <- rcs(log(simdata.tt$eventtime), knots=prep$iknots)
  simdata.tt$x1t1 <- simdata.tt$x1 * time.spline[ ,1]
  simdata.tt$x1t2 <- simdata.tt$x1 * time.spline[ ,2]
  simdata.tt$x2t1 <- simdata.tt$x2 * time.spline[ ,1]
  simdata.tt$x2t2 <- simdata.tt$x2 * time.spline[ ,2]
  simdata.tt$x3t1 <- simdata.tt$x3 * time.spline[ ,1]
  simdata.tt$x3t2 <- simdata.tt$x3 * time.spline[ ,2]
  simdata.tt$x4t1 <- simdata.tt$x4 * time.spline[ ,1]
  simdata.tt$x4t2 <- simdata.tt$x4 * time.spline[ ,2]
  simdata.tt$x5t1 <- simdata.tt$x5 * time.spline[ ,1]
  simdata.tt$x5t2 <- simdata.tt$x5 * time.spline[ ,2]
  simdata.tt$x6t1 <- simdata.tt$x6 * time.spline[ ,1]
  simdata.tt$x6t2 <- simdata.tt$x6 * time.spline[ ,2]
  simdata.tt$x7t1 <- simdata.tt$x7 * time.spline[ ,1]
  simdata.tt$x7t2 <- simdata.tt$x7 * time.spline[ ,2]
  simdata.tt$x8t1 <- simdata.tt$x8 * time.spline[ ,1]
  simdata.tt$x8t2 <- simdata.tt$x8 * time.spline[ ,2]
  simdata.tt$x9t1 <- simdata.tt$x9 * time.spline[ ,1]
  simdata.tt$x9t2 <- simdata.tt$x9 * time.spline[ ,2]
  simdata.tt$x10t1 <- simdata.tt$x10 * time.spline[ ,1]
  simdata.tt$x10t2 <- simdata.tt$x10 * time.spline[ ,2]
  rm(time.spline)
  
  if(!higherp){
    form.tt.alt <- formula(paste("Surv(tstart,eventtime, status) ~ ", paste0("x", 1:nvarobs, collapse="+"),
                                 "+", paste0("x",1:nvarobs, "t1", collapse="+"), "+", paste0("x",1:nvarobs, "t2", collapse="+")))
  } else {
    form.tt.alt <- formula(paste("Surv(tstart,eventtime, status) ~ ", paste0("x", 1:nvarobs, collapse="+"),
                                 "+", paste0("x",1:nvarobs, "t1", collapse="+"), "+", paste0("x",1:nvarobs, "t2", collapse="+"),
                                 "+", paste0("n",1:20, collapse="+")))
  }
  
  tt.alt <- coxph(form.tt.alt, data=simdata.tt, control=coxph.control(timefix=FALSE))

  time.per.model[s,4] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    val <- validation
    val$eventtime <- time.grid[ ,x]
    exp(-predict.cox.td(mod=tt.alt, newdata=val,
                        ndtimes=NULL,
                        type="expected",
                        iknots=prep$iknots,
                        logtime = TRUE))})

  rmse[s, 4, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 4, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[4]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[4]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[4]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[4]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        tt.alt, type="coxtt", logtime=TRUE, iknots=prep$iknots)$Ctd
  }))
  cstats.fixed[[s]][[4]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat, tt.alt)

  ###################################
  ## cox lasso                     ##
  ###################################

  model.time <- proc.time()

  cv <- cv.glmnet(x=as.matrix(simdata[ ,c(grep("x", names(simdata))[1:nvarobs],
                                          grep("n[0-9]", names(simdata)))]),
                  y=Surv(simdata$eventtime, simdata$status),
                  family="cox",
                  alpha=1,
                  nfolds=nfolds)

  cox.lasso <- coxph(form, data=simdata,
                     init=as.numeric(coef(cv, s=cv$lambda.min)[all.vars(form)[-c(1,2)], ]),
                     control=coxph.control(iter.max=0, timefix=FALSE))

  time.per.model[s,5] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    exp(-predict(cox.lasso, newdata=pred.data.other[[x]], type="expected"))
  })

  rmse[s, 5, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 5, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[5]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[5]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[5]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[5]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        cox.lasso, type="coxph")$Ctd
  }))
  cstats.fixed[[s]][[5]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat)

  ####################################
  ## cox ridge tt                   ##
  ####################################
  ## Possible since glmnet v4.1:
  ## https://statisticaloddsandends.wordpress.com/2021/01/14/glmnet-v4-1-regularized-cox-models-for-start-stop-and-stratified-data/

  model.time <- proc.time()
  
  # ptm <- proc.time()
  cv <- cv.glmnet(x=as.matrix(simdata.tt[ ,c(grep("x", names(simdata.tt))[1:nvarobs],
                                             grep("t1", names(simdata.tt)),
                                             grep("t2", names(simdata.tt)),
                                             grep("n[0-9]", names(simdata.tt)))]),
                  y=Surv(simdata.tt$tstart, simdata.tt$eventtime, simdata.tt$status),
                  family="cox",
                  alpha=0,
                  nfolds=nfolds)
  # proc.time() - ptm

  form.tt.ridge <- formula(paste0("Surv(tstart, eventtime, status) ~ ",
                                  paste0(names(simdata.tt)[c(grep("x", names(simdata.tt))[1:nvarobs],
                                                             grep("t1", names(simdata.tt)),
                                                             grep("t2", names(simdata.tt)),
                                                             grep("n[0-9]", names(simdata.tt)))], collapse="+")))

  tt.alt.ridge <- coxph(form.tt.ridge,
                        data=simdata.tt,
                        init=as.numeric(coef(cv, s=cv$lambda.min)[all.vars(form.tt.ridge)[-c(1:3)], ]),
                        control=coxph.control(iter.max=0, timefix=FALSE))

  time.per.model[s,6] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    val <- validation
    val$eventtime <- time.grid[ ,x]
    exp(-predict.cox.td(mod=tt.alt.ridge, newdata=val,
                        ndtimes=NULL,
                        type="expected",
                        logtime = TRUE,
                        iknots=prep$iknots))})

  rmse[s, 6, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 6, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[6]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[6]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[6]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[6]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        tt.alt.ridge, type="coxtt", logtime=TRUE, iknots=prep$iknots)$Ctd
  }))
  cstats.fixed[[s]][[6]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat, tt.alt.ridge)

  ####################################
  ## (p)rstpm2                      ##
  ####################################
  ## NB ## (p)rstpm2 models are on the logH scale and by default use log(time) (unaltered here)

  # fixed
  model.time <- proc.time()

  mod.stpm.ph <- stpm2(form,
                       smooth.formula=~ns(log(eventtime), 2),
                       data=simdata)

  time.per.model[s,7] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    predict(mod.stpm.ph, newdata=pred.data.other[[x]], type="surv")
  })
  rmse[s, 7, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 7, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[7]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[7]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[7]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[7]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        mod.stpm.ph, type="rstpm")$Ctd
  }))
  cstats.fixed[[s]][[7]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  # fixed nonPH
  model.time <- proc.time()

  smooth.formula <- formula(paste("~ rcs(log(eventtime), prep$knots) + ",
                                  paste0(all.vars(form)[grep("x", all.vars(form))], ":rcs(log(eventtime), prep$iknots)", collapse="+")))
  
  mod.stpm <- stpm2(form,
                    smooth.formula=smooth.formula,
                    data=simdata)

  time.per.model[s,8] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    predict(mod.stpm, newdata=pred.data.other[[x]], type="surv")
  })
  rmse[s, 8, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 8, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[8]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[8]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[8]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[8]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        mod.stpm.ph, type="rstpm")$Ctd
  }))
  cstats.fixed[[s]][[8]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat)

  # penalized (functional form)
  model.time <- proc.time()

  mod.pstpm.ph <- pstpm2(form,
                         smooth.formula= ~ s(log(eventtime), k=-1),
                         data=simdata)

  time.per.model[s,9] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    predict(mod.pstpm.ph, newdata=pred.data.other[[x]], type="surv")
  })
  rmse[s, 9, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 9, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[9]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[9]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[9]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[9]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        mod.stpm.ph, type="rstpm")$Ctd
  }))
  cstats.fixed[[s]][[9]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat)

  # penalized (functional form) nonPH
  model.time <- proc.time()

  smooth.formula <- formula(paste("~ s(log(eventtime), k=-1) + ", paste0("s(log(eventtime), by=", paste0("x",1:nvarobs), ", k=5)", collapse=" + ")))
  
  if(!higherp){
    form <- formula(Surv(eventtime, status) ~ 1)
  } else {
    form <- formula(paste("Surv(eventtime, status) ~ 1 +", paste(names(simdata)[grep("n[0-9]", names(simdata))], collapse="+"))) 
  }
  
  mod.pstpm <- pstpm2(form,
                      smooth.formula= smooth.formula,
                      data=simdata)

  time.per.model[s,10] <- (proc.time() - model.time)["elapsed"]

  Shat <- sapply(1:7, function(x){
    predict(mod.pstpm, newdata=pred.data.other[[x]], type="surv")
  })
  rmse[s, 10, ] <- apply(Shat - Strue, 2, function(x) sqrt(mean(x^2)))
  e90[s, 10, ] <- apply(Shat - Strue, 2, function(x) quantile(abs((x)), probs=.9, na.rm=TRUE))
  cal[[s]][[10]] <- tryCatch(
    expr=lapply(1:length(fixed.times), function(t)
      calibrate(Shat=Shat[ ,1+t],
                time=fixed.times[t], validation=validation,
                predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
    error=function(e) NULL)
  cstats[[s]][[10]] <- sapply(1:ncol(time.grid), function(t) rcorr.cens(Shat[ ,t], Strue[ ,t]))
  brier[[s]][[10]] <- bsctd(validation$eventtime, validation$status, fixed.times, Shat[ ,-1], weights = "msurv")
  cstatstd[[s]][[10]] <- mean(sapply(1:10, function(x){
    ss <- sample(1:nrow(validation), 100)
    ctd(validation$eventtime[ss], validation$status[ss],
        X=pred.data.other[[1]][ss ,-which(colnames(pred.data.other[[1]]) == "eventtime")],
        mod.stpm.ph, type="rstpm")$Ctd
  }))
  cstats.fixed[[s]][[10]] <- ctd.fixed(validation$eventtime, validation$status, fixed.times, Shat[ ,-1])

  rm(Shat)

  time.elapsed.sec[s] <- as.numeric(round((proc.time() - ptm)["elapsed"])); cleanup(ws.names)

  print(paste("Iteration",iter,"; Sample size setting", s))
}

results <- list(time.elapsed.sec=time.elapsed.sec,
                time.per.model.sec=time.per.model,
                rmse=rmse,
                e90=e90,
                cal=cal,
                predict.grid=predict.grid,
                cstats=cstats,
                brier=brier,
                cstatstd=cstatstd,
                cstats.fixed=cstats.fixed,
                coxzph.results=coxzph.results,
                time.per.model=time.per.model,
                sessionInfo=sessionInfo())
