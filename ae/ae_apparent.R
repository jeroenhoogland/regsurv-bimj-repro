source("sim/helperfunctions.R")
source("ae/helperfunctions_ae.R")

time.per.model <- rep(NA, 10)

# prepare data for easier predictions lateron
pred.data.regsurv <- list(as.matrix(data.frame("eventtime"=veteran$eventtime, X)))
pred.data.regsurv[2:ncol(time.grid)] <- lapply(1:length(fixed.times), function(t){
  as.matrix(data.frame("eventtime"=
                         rep(fixed.times[t], 
                             length(veteran$eventtime)), X))})

pred.data.other <- list(veteran)
pred.data.other[2:ncol(time.grid)] <- lapply(1:length(fixed.times), function(t){
  val.fixed <- veteran
  val.fixed$eventtime <- fixed.times[t]
  val.fixed})

predict.grid <- cbind(seq(.2,.9,length.out=51),
                      seq(.3,1,length.out=51),
                      seq(.4,1,length.out=51),
                      seq(.45,1,length.out=51),
                      seq(.5,1,length.out=51))
predict.grid <- predict.grid[-51, ]

nfolds <- 10
maxit <- 500L
brier <- cstats <- cstats.fixed <- cal <- vector(mode = "list", length = 12)
coxzph.results <- list()

ntimebasis <- 4
nitimebasis <- 1


###################################
## log time log hazard           ##
###################################

model.time <- proc.time()

prep <- survprep(tte=as.numeric(veteran$eventtime),
                 delta=as.numeric(veteran$status),
                 X=as.matrix(X),
                 model.scale="loghazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=1:ncol(X),
                 nitimebasis=nitimebasis,
                 qpoints=30)

groups <- list()
groups[[1]] <- 1
groups[[2]] <- prep$which.param[[1]][-1]
groups[[3]]  <- prep$which.param[[2]][1]
groups[[4]]  <- prep$which.param[[2]][2:4]
groups[[5]]  <- prep$which.param[[2]][5]
groups[[6]]  <- prep$which.param[[2]][6]
groups[[7]]  <- prep$which.param[[2]][7]
groups[[8]]  <- prep$which.param[[2]][8]
groups[[9]] <- prep$which.param[[3]][1] # nonPH trt
groups[[10]] <- prep$which.param[[3]][2:4] # nonPH celltype
groups[[11]] <- prep$which.param[[3]][5] # nonPH karno
groups[[12]] <- prep$which.param[[3]][6] # nonPH diagtime
groups[[13]] <- prep$which.param[[3]][7] # nonPH age
groups[[14]] <- prep$which.param[[3]][8] # nonPH prior

penpars <- c(0, rep(1, length(groups)-1)) #  penalize all except intercept
l1l2 <- c(rep(0,2),rep(1, length(groups)-2)) # ridge baseline, (group) lasso for the rest

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

time.per.model[1] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                      newdata=pred.data.regsurv[[x]], type="surv"))

# Brier score td
brier[[1]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[1]] <- ctd(veteran$eventtime, veteran$status, X=X, obj$mod, type="regsurv",
                   prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed[[1]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status, X=X, obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal[[1]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat)
rm(obj)

###################################
## log time log Hazard           ##
###################################

model.time <- proc.time()

prep <- survprep(tte=as.numeric(veteran$eventtime),
                 delta=as.numeric(veteran$status),
                 X=as.matrix(X),
                 model.scale="logHazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=1:ncol(X),
                 nitimebasis=nitimebasis,
                 qpoints=NULL)

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

time.per.model[2] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                      newdata=pred.data.regsurv[[x]], type="surv"))

# Brier score td
brier[[2]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[2]] <- ctd(veteran$eventtime, veteran$status, X=X, obj$mod, type="regsurv",
                   prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed[[2]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status, X=X, obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal[[2]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat)
rm(obj)

###################################
## coxph                         ##
###################################

model.time <- proc.time()

form <- formula(Surv(eventtime, status) ~ trt + celltype + karno +
                  diagtime + age + prior)

# see coxph.control() help for timefix (especially in simulated data without ties)
mod.coxph <- coxph(form, data=veteran, control=coxph.control(timefix=FALSE))

time.per.model[3] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  exp(-predict(mod.coxph, newdata=pred.data.other[[x]], type="expected"))
})

# Brier score td
brier[[3]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[3]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   mod.coxph, type="coxph")

cstats.fixed[[3]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      mod.coxph, type="coxph")
})

# Calibration
cal[[3]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, mod.coxph)

###################################
## coxph.tt                      ##
###################################

model.time <- proc.time()

# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
dtimes <- unique(veteran$eventtime[veteran$status==1])
data.tt <- survSplit(Surv(eventtime, status==1) ~., veteran, cut=dtimes)
names(data.tt)[names(data.tt) == "event"] <- "status"
data.tt$time.spline <- log(data.tt$eventtime)
mmtt <- model.matrix( ~  tstart + eventtime + status + trt + celltype + karno +
                        diagtime + age + prior + (trt + celltype + karno +
                                                    diagtime + age + prior):time.spline,
                      data=data.tt)
mmtt <- mmtt[ ,-1]
vars <- c("trt", "celltypesmallcell", "celltypeadeno", "celltypelarge", "karno", "diagtime", "age", "prior")
mmtt <- data.frame(mmtt)
names(mmtt) <- c( "tstart", "eventtime", "status",vars,paste0(vars, "t1"))
form.tt.alt <- formula(paste("Surv(tstart,eventtime, status) ~ ", paste0(names(mmtt)[-c(1:3)], collapse = "+")))
tt.alt <- coxph(form.tt.alt, data=mmtt, control=coxph.control(timefix=FALSE))

time.per.model[4] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  val <- veteran
  val$eventtime <- time.grid[ ,x]
  exp(-predict.cox.td(mod=tt.alt, newdata=val,
                      ndtimes=NULL,
                      type="expected",
                      logtime = TRUE))})

# Brier score td
brier[[4]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[4]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   tt.alt, type="coxtt", logtime = TRUE, iknots = NULL)

cstats.fixed[[4]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      tt.alt, type="coxtt", logtime = TRUE, iknots = NULL)
})

# Calibration
cal[[4]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, tt.alt)

###################################
## cox ridge                     ##
## (due to grouped celltype)     ##
###################################

model.time <- proc.time()

cv <- cv.glmnet(x=X,
                y=Surv(veteran$eventtime, veteran$status),
                family="cox",
                alpha=0,
                nfolds=nfolds)

cox.ridge <- coxph(form, data=veteran,
                   init=as.numeric(coef(cv, s=cv$lambda.min)),
                   control=coxph.control(iter.max=0, timefix=FALSE))

time.per.model[5] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  exp(-predict(cox.ridge, newdata=pred.data.other[[x]], type="expected"))
})

# Brier score td
brier[[5]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[5]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   cox.ridge, type="coxph")

cstats.fixed[[5]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      cox.ridge, type="coxph")
})

# Calibration
cal[[5]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, cox.ridge)

######################################
## cox ridge tt (does not converge) ##
######################################

# # https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
# # dtimes <- unique(veteran$eventtime[veteran$status==1])
# # dtimes <- unique(quantile(veteran$eventtime[veteran$status==1], probs=seq(.05,.95, .05), type=1))
# # dtimes <- unique(quantile(veteran$eventtime[veteran$status==1], probs=seq(.1,.9, .1), type=1))
# dtimes <- unique(quantile(veteran$eventtime[veteran$status==1], probs=.5, type=1))
# data.tt <- survSplit(Surv(eventtime, status==1) ~., veteran, cut=dtimes)
# names(data.tt)[names(data.tt) == "event"] <- "status"
# data.tt$time.spline <- log(data.tt$eventtime)
# mmtt <- model.matrix( ~  tstart + eventtime + status + trt + celltype + karno +
#                       diagtime + age + prior + (trt + celltype + karno +
#                       diagtime + age + prior):time.spline,
#                       data=data.tt)
# mmtt <- mmtt[ ,-1]
# vars <- c("trt", "celltypesmallcell", "celltypeadeno", "celltypelarge", "karno", "diagtime", "age", "prior")
# mmtt <- data.frame(mmtt)
# names(mmtt) <- c( "tstart", "eventtime", "status",vars,paste0(vars, "t1"))
#
# cv <- cv.glmnet(x=as.matrix(mmtt[ ,-c(1:3)]),
#                 y=Surv(mmtt$tstart, mmtt$eventtime, mmtt$status),
#                 family="cox",
#                 alpha=0,
#                 nfolds=nfolds)
#
# form.tt.ridge <- formula(paste("Surv(tstart,eventtime, status) ~ ", paste0(names(mmtt)[-c(1:3)], collapse = "+")))
# tt.alt.ridge <- coxph(form.tt.ridge,
#                       data=mmtt,
#                       init=as.numeric(coef(cv, s=cv$lambda.min)),
#                       control=coxph.control(iter.max=0, timefix=FALSE))

####################################
## (p)rstpm2                      ##
####################################

model.time <- proc.time()

## NB ## (p)rstpm2 models are on the logH scale and by default use log(time) (unaltered here)
# fixed
mod.stpm.ph <- stpm2(form,
                     smooth.formula=~ns(log(eventtime), 2),
                     data=veteran)

time.per.model[7] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  predict(mod.stpm.ph, newdata=pred.data.other[[x]], type="surv")
})

# Brier score td
brier[[7]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[7]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   mod.stpm.ph, type="rstpm")

cstats.fixed[[7]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      mod.stpm.ph, type="rstpm")
})

# Calibration
cal[[7]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, mod.stpm.ph)

# fixed nonPH
model.time <- proc.time()
smooth.formula <- formula(paste("~ log(eventtime) + ",
                                paste0(all.vars(form)[-c(1,2)], ":log(eventtime)", collapse="+")))
mod.stpm <- stpm2(form,
                  smooth.formula=smooth.formula,
                  data=veteran)

time.per.model[8] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  predict(mod.stpm, newdata=pred.data.other[[x]], type="surv")
})

# Brier score td
brier[[8]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[8]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   mod.stpm, type="rstpm")

cstats.fixed[[8]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      mod.stpm, type="rstpm")
})

# Calibration
cal[[8]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, mod.stpm)

# penalized (functional form)
model.time <- proc.time()
mod.pstpm.ph <- pstpm2(form,
                       smooth.formula= ~ s(log(eventtime), k=-1),
                       data=veteran)

time.per.model[9] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  predict(mod.pstpm.ph, newdata=pred.data.other[[x]], type="surv")
})

# Brier score td
brier[[9]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[9]] <- ctd(veteran$eventtime, veteran$status,
                   X=veteran[ ,-which(names(veteran)=="eventtime")],
                   mod.pstpm.ph, type="rstpm")

cstats.fixed[[9]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      mod.pstpm.ph, type="rstpm")
})

# Calibration
cal[[9]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, mod.pstpm.ph)

# penalized (functional form) nonPH
model.time <- proc.time()
mod.pstpm <- pstpm2(Surv(eventtime, status) ~ trt + celltype +
                      karno + diagtime + age + prior,
                    smooth.formula= ~s(log(eventtime), k=-1) + trt:log(eventtime) +
                      celltype:log(eventtime) + karno:log(eventtime) + diagtime:log(eventtime) +
                      age:log(eventtime) + prior:log(eventtime),
                    data=veteran)

time.per.model[10] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x){
  predict(mod.pstpm, newdata=pred.data.other[[x]], type="surv")
})

# Brier score td
brier[[10]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[10]] <- ctd(veteran$eventtime, veteran$status,
                    X=veteran[ ,-which(names(veteran)=="eventtime")],
                    mod.pstpm, type="rstpm")

cstats.fixed[[10]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status,
      X=veteran[ ,-which(names(veteran)=="eventtime")],
      mod.pstpm, type="rstpm")
})

# Calibration
cal[[10]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat, mod.pstpm)

###################################
## log time log hazard           ##
###################################

model.time <- proc.time()
prep <- survprep(tte=as.numeric(veteran$eventtime),
                 delta=as.numeric(veteran$status),
                 X=as.matrix(X),
                 model.scale="loghazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=NULL,
                 nitimebasis=NULL,
                 qpoints=30)

groups <- list()
groups[[1]] <- 1
groups[[2]] <- prep$which.param[[1]][-1]
groups[[3]]  <- prep$which.param[[2]][1]
groups[[4]]  <- prep$which.param[[2]][2:4]
groups[[5]]  <- prep$which.param[[2]][5]
groups[[6]]  <- prep$which.param[[2]][6]
groups[[7]]  <- prep$which.param[[2]][7]
groups[[8]]  <- prep$which.param[[2]][8]

penpars <- c(0,0, rep(1, length(groups)-2)) #  unpenalized baseline
l1l2 <- c(rep(0,2),rep(1, length(groups)-2)) # ridge baseline, (group) lasso for the rest

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

time.per.model[11] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                      newdata=pred.data.regsurv[[x]], type="surv"))

# Brier score td
brier[[11]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[11]] <- ctd(veteran$eventtime, veteran$status, X=X, obj$mod, type="regsurv",
                    prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed[[11]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status, X=X, obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal[[11]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat)
rm(obj)

###################################
## log time log Hazard           ##
###################################

model.time <- proc.time()

prep <- survprep(tte=as.numeric(veteran$eventtime),
                 delta=as.numeric(veteran$status),
                 X=as.matrix(X),
                 model.scale="loghazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=NULL,
                 nitimebasis=NULL,
                 qpoints=30)

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

time.per.model[12] <- (proc.time() - model.time)["elapsed"]

Shat <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                      newdata=pred.data.regsurv[[x]], type="surv"))

# Brier score td
brier[[12]] <- bsctd(veteran$eventtime, veteran$status, fixed.times, Shat[ ,-1], weights = "msurv")

# Ctd
cstats[[12]] <- ctd(veteran$eventtime, veteran$status, X=X, obj$mod, type="regsurv",
                    prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed[[12]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime, v.censored$status, X=X, obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal[[12]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat[ ,1+t],
              time=fixed.times[t], validation=veteran,
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=TRUE)),
  error=function(e) NULL)

rm(Shat)
rm(obj)

apparent <- list(brier=brier,
                 cstats=cstats,
                 cstats.fixed=cstats.fixed,
                 cal=cal)
