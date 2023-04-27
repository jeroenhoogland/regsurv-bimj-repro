# laod packages and data
library(regsurv)
library(survival)
library(polspline)
library(rstpm2)
library(glmnet)
library(Hmisc)
library(dplyr)

source("sim/helperFunctions.R") 
source("ae/helperfunctions_ae.R")

ptm <- proc.time()

index <- sample(1:nrow(veteran), size = nrow(veteran), replace=TRUE)
out <- (1:nrow(veteran))[-index]

names(veteran)[which(names(veteran) == "time")] <- "eventtime"

train <- veteran[index, ]
test <- veteran[out, ]

# model matrix
mm <- model.matrix( ~ trt + celltype + karno +
                      diagtime + age + prior,
                    data=veteran)
X <- mm[ ,-1]
fixed.times <- c(seq(60, 180, 30))
time.grid <- cbind(veteran$eventtime, do.call(cbind, lapply(fixed.times, function(x) rep(x, nrow(veteran)))))

# some preparation for predictions
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

nfolds <- 10 # for parameter tuning
maxit <- 500L # for regsurv
brier.oob <- cstats.oob <- cstats.fixed.oob <- cal.oob <- vector(mode = "list", length = 12)

ntimebasis <- 4
nitimebasis <- 1

###################################
## log time log hazard           ##
###################################

prep <- survprep(tte=as.numeric(train$eventtime),
                 delta=as.numeric(train$status),
                 X=as.matrix(X[index, ]),
                 model.scale="loghazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=1:ncol(X[index, ]),
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
l1l2 <- c(rep(0,length(groups)-6),rep(1, 6)) # as in sims

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

Shat.oob <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                          newdata=pred.data.regsurv[[x]][-index, ], type="surv"))

# Brier score td
brier.oob[[1]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[1]] <- ctd(veteran$eventtime[-index], veteran$status[-index], X=X[-index, ], obj$mod, type="regsurv",
                       prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed.oob[[1]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index], X=X[-index, ], obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal.oob[[1]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob)
rm(obj)

###################################
## log time log Hazard           ##
###################################

prep <- survprep(tte=as.numeric(train$eventtime),
                 delta=as.numeric(train$status),
                 X=as.matrix(X[index, ]),
                 model.scale="logHazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=1:ncol(X),
                 nitimebasis=nitimebasis,
                 qpoints=NULL)

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

Shat.oob <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                          newdata=pred.data.regsurv[[x]][-index, ], type="surv"))

# Brier score td
brier.oob[[2]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[2]] <- ctd(veteran$eventtime[-index], veteran$status[-index], X=X[-index, ], obj$mod, type="regsurv",
                       prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed.oob[[2]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index], X=X[-index, ], obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal.oob[[2]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob)
rm(obj)


###################################
## coxph                         ##
###################################

form <- formula(Surv(eventtime, status) ~ trt + celltype + karno +
                  diagtime + age + prior)

# see coxph.control() help for timefix (especially in simulated data without ties)
mod.coxph <- coxph(form, data=train, control=coxph.control(timefix=FALSE))

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  exp(-predict(mod.coxph, newdata=pred.data.other[[x]][-index, ], type="expected"))
})

# Brier score td
brier.oob[[3]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[3]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       mod.coxph, type="coxph")

cstats.fixed.oob[[3]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      mod.coxph, type="coxph")
})

# Calibration
cal.oob[[3]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, mod.coxph)

###################################
## coxph.tt                      ##
###################################

# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
dtimes <- unique(train$eventtime[train$status==1])
data.tt <- survSplit(Surv(eventtime, status==1) ~., train, cut=dtimes)
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

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  val <- veteran[-index, ]
  val$eventtime <- time.grid[-index ,x]
  exp(-predict.cox.td(mod=tt.alt, newdata=val,
                      ndtimes=NULL,
                      type="expected",
                      logtime = TRUE))})

# Brier score td
brier.oob[[4]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[4]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       tt.alt, type="coxtt", logtime=TRUE)

cstats.fixed.oob[[4]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      tt.alt, type="coxtt", logtime=TRUE)
})

# Calibration
cal.oob[[4]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, tt.alt)

###################################
## cox ridge                     ##
## (due to grouped celltype)     ##
###################################

cv <- cv.glmnet(x=X[index, ],
                y=Surv(train$eventtime, train$status),
                family="cox",
                alpha=0,
                nfolds=nfolds)

cox.ridge <- coxph(form, data=train,
                   init=as.numeric(coef(cv, s=cv$lambda.min)),
                   control=coxph.control(iter.max=0, timefix=FALSE))

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  exp(-predict(cox.ridge, newdata=pred.data.other[[x]][-index, ], type="expected"))
})

# Brier score td
brier.oob[[5]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[5]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       cox.ridge, type="coxph")

cstats.fixed.oob[[5]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      cox.ridge, type="coxph")
})

# Calibration
cal.oob[[5]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, cox.ridge)

######################################
## cox ridge tt (does not converge) ##
######################################

# # https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
# # dtimes <- unique(veteran$eventtime[veteran$status==1])
# dtimes <- unique(quantile(veteran$eventtime[veteran$status==1], probs=seq(0,1,length.out=11), type=1))
# data.tt <- survSplit(Surv(eventtime, status==1) ~., veteran, cut=dtimes)
# names(data.tt)[names(data.tt) == "event"] <- "status"
# data.tt$time.spline <- data.tt$eventtime
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
## NB ## (p)rstpm2 models are on the logH scale and by default use log(time) (unaltered here)
# fixed
mod.stpm.ph <- stpm2(form,
                     smooth.formula=~ns(log(eventtime), 2),
                     data=train)

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  predict(mod.stpm.ph, newdata=pred.data.other[[x]][-index, ], type="surv")
})

# Brier score td
brier.oob[[7]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[7]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       mod.stpm.ph, type="rstpm")

cstats.fixed.oob[[7]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      mod.stpm.ph, type="rstpm")
})

# Calibration
cal.oob[[7]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, mod.stpm.ph)

# fixed nonPH
smooth.formula <- formula(paste("~ log(eventtime) + ",
                                paste0(all.vars(form)[-c(1,2)], ":log(eventtime)", collapse="+")))
mod.stpm <- stpm2(form,
                  smooth.formula=smooth.formula,
                  data=train)

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  predict(mod.stpm, newdata=pred.data.other[[x]][-index, ], type="surv")
})

# Brier score td
brier.oob[[8]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[8]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       mod.stpm, type="rstpm")

cstats.fixed.oob[[8]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      mod.stpm, type="rstpm")
})

# Calibration
cal.oob[[8]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, mod.stpm)

# penalized (functional form)
mod.pstpm.ph <- pstpm2(form,
                       smooth.formula= ~ s(log(eventtime), k=-1),
                       data=train)

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  predict(mod.pstpm.ph, newdata=pred.data.other[[x]][-index, ], type="surv")
})

# Brier score td
brier.oob[[9]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[9]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                       X=veteran[-index ,-which(names(veteran)=="eventtime")],
                       mod.pstpm.ph, type="rstpm")

cstats.fixed.oob[[9]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      mod.pstpm.ph, type="rstpm")
})

# Calibration
cal.oob[[9]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob, mod.pstpm.ph)

# penalized (functional form) nonPH
mod.pstpm <- pstpm2(Surv(eventtime, status) ~ trt + celltype +
                      karno + diagtime + age + prior,
                    smooth.formula= ~s(log(eventtime), k=-1) + trt:log(eventtime) +
                      celltype:log(eventtime) + karno:log(eventtime) + diagtime:log(eventtime) +
                      age:log(eventtime) + prior:log(eventtime),
                    data=train)

Shat.oob <- sapply(1:ncol(time.grid), function(x){
  predict(mod.pstpm, newdata=pred.data.other[[x]][-index, ], type="surv")
})

# Brier score td
brier.oob[[10]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[10]] <- ctd(veteran$eventtime[-index], veteran$status[-index],
                        X=veteran[-index ,-which(names(veteran)=="eventtime")],
                        mod.pstpm, type="rstpm")

cstats.fixed.oob[[10]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index],
      X=veteran[-index ,-which(names(veteran)=="eventtime")],
      mod.pstpm, type="rstpm")
})

# Calibration
cal.oob[[10]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)


rm(Shat.oob, mod.pstpm)

## log time log hazard PH ##
prep <- survprep(tte=as.numeric(train$eventtime),
                 delta=as.numeric(train$status),
                 X=as.matrix(X[index, ]),
                 model.scale="loghazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=NULL,
                 nitimebasis=NULL,
                 qpoints=30)

# prep$which.param
# head(prep$mm$d)
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

Shat.oob <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                          newdata=pred.data.regsurv[[x]][-index, ], type="surv"))

# Brier score td
brier.oob[[11]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[11]] <- ctd(veteran$eventtime[-index], veteran$status[-index], X=X[-index, ], obj$mod, type="regsurv",
                        prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed.oob[[11]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index], X=X[-index, ], obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal.oob[[11]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob)
rm(obj)

## log time log Hazard PH  ##
prep <- survprep(tte=as.numeric(train$eventtime),
                 delta=as.numeric(train$status),
                 X=as.matrix(X[index, ]),
                 model.scale="logHazard",
                 time.scale="logtime",
                 spline.type="rcs",
                 ntimebasis=ntimebasis,
                 tv=NULL,
                 nitimebasis=NULL,
                 qpoints=30)

obj <- group.rs.pred(prep=prep, penpars=penpars, l1l2=l1l2, groups=groups, nfolds=nfolds, maxit=maxit)

Shat.oob <- sapply(1:ncol(time.grid), function(x) predict(obj$mod, prep, lambda.index=obj$cv$lambda.min.index,
                                                          newdata=pred.data.regsurv[[x]][-index, ], type="surv"))

# Brier score td
brier.oob[[12]] <- bsctd(veteran$eventtime[-index], veteran$status[-index], fixed.times, Shat.oob[ ,-1], weights = "msurv")

# Ctd
cstats.oob[[12]] <- ctd(veteran$eventtime[-index], veteran$status[-index], X=X[-index, ], obj$mod, type="regsurv",
                        prep=prep, lambda.index=obj$cv$lambda.min.index)

cstats.fixed.oob[[12]] <- sapply(fixed.times, function(x){
  v.censored <- veteran
  v.censored$status <- ifelse(v.censored$eventtime > x, 0, v.censored$status)
  v.censored$eventtime <- ifelse(v.censored$eventtime > x, x, v.censored$eventtime)
  ctd(v.censored$eventtime[-index], v.censored$status[-index], X=X[-index, ], obj$mod, type="regsurv",
      prep=prep, lambda.index=obj$cv$lambda.min.index)
})

# Calibration
cal.oob[[12]] <-  tryCatch(
  expr=lapply(1:length(fixed.times), function(t)
    calibrate(Shat=Shat.oob[ ,1+t],
              time=as.numeric(fixed.times[t]), validation=veteran[-index, ],
              predict.grid=predict.grid[ ,t], plot=FALSE, store.predict.grid=FALSE)),
  error=function(e) NULL)

rm(Shat.oob)
rm(obj)

time.elapsed.sec <- (proc.time() - ptm)["elapsed"]

boot <- list(time.elapsed.sec=time.elapsed.sec,
             brier.oob=brier.oob,
             cstats.oob=cstats.oob,
             cstats.fixed.oob=cstats.fixed.oob,
             cal.oob=cal.oob)