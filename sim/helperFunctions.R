# Helper functions for the simulation R script

# regsurv with tryCatch for simulation runs
group.rs.pred <- function(prep, penpars, l1l2, groups, nfolds, lambda.grid=NULL,
                          maxit=maxit, pred.data, print=FALSE){
  mod <- tryCatch(
    expr=regsurv(prep, penpars, l1l2, groups, lambda.grid=lambda.grid,
                 maxit=maxit, print=print),
    error=function(e) NULL)
  
  if(!is.null(mod)){
    cv <- tryCatch(
      expr=cv.regsurv(mod, prep, plot=FALSE, print=print, nfolds=nfolds),
      error=function(e) NULL)
  } else {
    cv <- NULL
  }
  
  if(!is.null(mod) & !is.null(cv)){
    return(list(mod=mod, cv=cv))
  } else {
    return(NULL)
  }
}

# Calibration according to Austin 2020
# validation data should have cols eventtime and status
calibrate <- function(Shat, time, validation, predict.grid=NULL, plot=FALSE, store.predict.grid=FALSE){
  
  # coxph calibration model
  Fhat <- 1-Shat
  Fhat.cll <- log(-log(1-Fhat))
  
  # restrict Fhat.cll to finite range
  Fhat.cll <- ifelse(Fhat.cll == Inf, max(Fhat.cll[is.finite(Fhat.cll)]), Fhat.cll)
  Fhat.cll <- ifelse(Fhat.cll == -Inf, min(Fhat.cll[is.finite(Fhat.cll)]), Fhat.cll)
  
  kn <- stats::quantile(Fhat.cll, probs = c(.1, 0.5, .9)) # as suggested by Austin 2020
  
  calibrate.cox <- tryCatch(
    expr=coxph(Surv(eventtime, status) ~ rcs(Fhat.cll, knots = kn), x=TRUE, data=validation),
    error=function(e) NULL)
  
  if(is.null(predict.grid)){
    predict.grid <- seq(quantile(Fhat,probs=0.01),
                        quantile(Fhat,probs=0.99),length=100)
  }
  predict.grid.cll <- log(-log(1-predict.grid))
  predict.grid.df <- data.frame(predict.grid)
  predict.grid.cll.df <- data.frame(predict.grid.cll)
  names(predict.grid.df) <- "Fhat"
  names(predict.grid.cll.df) <- "Fhat.cll"
  
  if(!is.null(calibrate.cox)){
    predict.calibrate.cox <- 1 - exp(-predict(calibrate.cox,
                                              newdata=data.frame(eventtime=time,status=1,Fhat.cll=predict.grid.cll.df),
                                              type="expected"))
    if(plot){
      plot(predict.grid,predict.calibrate.cox,type="l",lty=1,col="red",
           xlim=c(0,1),ylim=c(0,1),
           xlab = "Predicted probability of 1-year mortality",
           ylab = "Observed probability of 1-year mortality")
      par(new=T)
      plot(density(Fhat),axes=F,xlab=NA,ylab=NA,main="",xlim=c(0,1))
      axis(side=4)
    }
    
    cox.cal.stats <- list(ICI=mean(abs(predict.grid - predict.calibrate.cox), na.rm=TRUE),
                          E50=median(abs(predict.grid - predict.calibrate.cox), na.rm=TRUE),
                          E90=quantile(abs(predict.grid - predict.calibrate.cox), probs=.9, na.rm=TRUE),
                          cox.zph=tryCatch(
                            expr=cox.zph(calibrate.cox)$table,
                            error=function(e) NULL),
                          warn=NULL)
  } else {
    cox.cal.stats <- list(ICI=NULL,
                          E50=NULL,
                          E90=NULL,
                          cox.zph=NULL,
                          warn=c("(near) constant predictions with sd of Shat"=sd(Shat)))
    predict.calibrate.cox <- NULL
  }
  
  # hare calibration model
  calibrate.hare <- hare(data=validation$eventtime, delta=validation$status, cov=as.matrix(Fhat.cll))
  predict.calibrate.hare <- phare(time, predict.grid.cll.df, calibrate.hare)
  if(plot){
    plot(predict.grid,predict.calibrate.hare,type="l",lty=1,col="red",
         xlim=c(0,1),ylim=c(0,1),
         xlab = "Predicted probability of 1-year mortality",
         ylab = "Observed probability of 1-year mortality")
    par(new=T)
    plot(density(Fhat),axes=F,xlab=NA,ylab=NA,main="",xlim=c(0,1))
    axis(side=4)
  }
  
  hare.cal.stats <- list(ICI=mean(abs(predict.grid - predict.calibrate.hare), na.rm=TRUE),
                         E50=median(abs(predict.grid - predict.calibrate.hare), na.rm=TRUE),
                         E90=quantile(abs(predict.grid - predict.calibrate.hare), probs=.9, na.rm=TRUE))
  
  if(store.predict.grid){
    return(list(cox.cal.stats=cox.cal.stats,
                cox.cal.obs=predict.calibrate.cox,
                hare.cal.stats=hare.cal.stats,
                hare.cal.obs=predict.calibrate.hare,
                predict.grid=predict.grid))
  } else {
    return(list(cox.cal.stats=cox.cal.stats,
                cox.cal.obs=predict.calibrate.cox,
                hare.cal.stats=hare.cal.stats,
                hare.cal.obs=predict.calibrate.hare))
  }
}

# Helper function for time-dependent coefficients for the cox models in the same
# specification as for the other models.
# Currently not suited for left-truncation (not needed for current analyses)
predict.cox.td <- function(mod, newdata, ndtimes=NULL, type="expected", iknots, logtime=FALSE){
  if(is.null(ndtimes)){
    bh <- basehaz(mod, centered = FALSE)
    H0hat <- c(0, bh$hazard)
    Htimes <- c(0, bh$time)
  } else {
    bh <- basehaz(mod, centered = FALSE)
    dtimes <- unique(quantile(bh$time[-length(bh$time)], probs=seq(0,1,length.out=ndtimes), type=1))
    bh <- bh[bh$time %in% dtimes, ]
    H0hat <- c(0, bh$hazard)
    Htimes <- c(0, bh$time)
  }
  pred.matrix <- survSplit(Surv(eventtime, status==1) ~., newdata, cut=Htimes)
  names(pred.matrix)[names(pred.matrix) == "event"] <- "status"
  if(logtime){
    time.spline <- rcs(log(pred.matrix$eventtime), knots=iknots)
  } else {
    time.spline <- rcs(pred.matrix$eventtime, knots=iknots)
  }
  pred.matrix$x1t1 <- pred.matrix$x1 * time.spline[ ,1]
  pred.matrix$x1t2 <- pred.matrix$x1 * time.spline[ ,2]
  pred.matrix$x2t1 <- pred.matrix$x2 * time.spline[ ,1]
  pred.matrix$x2t2 <- pred.matrix$x2 * time.spline[ ,2]
  pred.matrix$x3t1 <- pred.matrix$x3 * time.spline[ ,1]
  pred.matrix$x3t2 <- pred.matrix$x3 * time.spline[ ,2]
  pred.matrix$x4t1 <- pred.matrix$x4 * time.spline[ ,1]
  pred.matrix$x4t2 <- pred.matrix$x4 * time.spline[ ,2]
  pred.matrix$x5t1 <- pred.matrix$x5 * time.spline[ ,1]
  pred.matrix$x5t2 <- pred.matrix$x5 * time.spline[ ,2]
  pred.matrix$x6t1 <- pred.matrix$x6 * time.spline[ ,1]
  pred.matrix$x6t2 <- pred.matrix$x6 * time.spline[ ,2]
  pred.matrix$x7t1 <- pred.matrix$x7 * time.spline[ ,1]
  pred.matrix$x7t2 <- pred.matrix$x7 * time.spline[ ,2]
  pred.matrix$x8t1 <- pred.matrix$x8 * time.spline[ ,1]
  pred.matrix$x8t2 <- pred.matrix$x8 * time.spline[ ,2]
  pred.matrix$x9t1 <- pred.matrix$x9 * time.spline[ ,1]
  pred.matrix$x9t2 <- pred.matrix$x9 * time.spline[ ,2]
  pred.matrix$x10t1 <- pred.matrix$x10 * time.spline[ ,1]
  pred.matrix$x10t2 <- pred.matrix$x10 * time.spline[ ,2]
  
  pred.matrix$H0hat <- predict(mod, newdata=pred.matrix, type="expected")
  H0hat <- pred.matrix %>% group_by(id) %>%
    summarise(cumhaz=sum(H0hat))
  H0hat <- as.numeric(H0hat$cumhaz)
  H0hat[is.na(H0hat)] <- 0 
  H0hat
}

# Time-dependent Brier score; independent censoring (according to Graf 1999)
bsctd <- function(eventtime, status, times, Shat, weights=c("linear", "msurv")){
  G <- coxph(Surv(eventtime, ifelse(status==1,0,1)) ~ 1)
  bh <- basehaz(G)
  br <- sapply(1:length(times), function(x){
    mean(Shat[ ,x]^2 * ifelse(eventtime <= times[x] & status == 1, 1, 0) *
           (1/exp(-predict(G, type="expected"))) +
           (1-Shat[ ,x])^2 * ifelse(eventtime > times[x], 1, 0) *
           (1/exp(-bh$hazard[max(which(bh$time <= times[x]))])))
  })
  
  if(weights == "linear"){
    int <- sapply(1:length(times), function(x){
      w <- times[1:x] / times[x]
      c(sum(br[1:x] * w), sum(br[1:x] * w) / sum(w))})
  } else if(weights == "msurv"){
    S <- coxph(Surv(eventtime, ifelse(status==1,1,0)) ~ 1)
    bh <- basehaz(S)
    St <- sapply(times, function(x) exp(-bh$hazard[max(which(bh$time <= x))]))
    int <- sapply(1:length(times), function(x){
      w <- (1-St[1:x]) / (1-exp(-bh$hazard[max(which(bh$time <= times[x]))]))
      c(sum(br[1:x] * w), sum(br[1:x] * w) / sum(w))})
  }
  data.frame("times"=times,
             "brctd"=br,
             "integrated"=int[1,],
             "averaged"=int[2,])
}

# Ctd (time-dependent C-statistic according to Antolini 2005)
ctd <- function(eventtime, status, X, model, type=c("regsurv", "coxph", "coxtt", "rstpm"),
                prep=NULL, lambda.index=NULL, logtime=NULL, iknots=NULL){
  cases <- which(status == 1)
  conc <- comp <- matrix(0, length(eventtime), length(eventtime))
  for(i in 1:length(cases)){
    comparable <- which(eventtime > eventtime[cases[i]] | (eventtime == eventtime[cases[i]] & status %in% 0))
    if(length(comparable) > 0){
      comp[i, comparable] <- 1
      if(type=="regsurv"){
        Stihat <- predict(model, prep, lambda.index,
                          newdata=as.matrix(data.frame("eventtime"=
                                                         rep(eventtime[cases[1]], length(comparable)+1),
                                                       X[c(cases[i],comparable), ])),
                          type="surv")
      }
      if(type=="coxph"){
        Stihat <- exp(-predict(model,
                               newdata=data.frame("eventtime"=rep(eventtime[cases[1]], length(comparable)+1),
                                                  X[c(cases[i],comparable), ]),
                               type="expected"))
      }
      if(type=="coxtt"){
        if(logtime){
          Stihat <- exp(-predict.cox.td(mod=model,
                                        newdata=data.frame("eventtime"=rep(eventtime[cases[1]], length(comparable)+1),
                                                           X[c(cases[i],comparable), ]),
                                        ndtimes=NULL, iknots = iknots,
                                        type="expected", logtime = TRUE))
        } else {
          Stihat <- exp(-predict.cox.td(mod=model,
                                        newdata=data.frame("eventtime"=rep(eventtime[cases[1]], length(comparable)+1),
                                                           X[c(cases[i],comparable), ]),
                                        ndtimes=NULL, iknots = iknots,
                                        type="expected"))
        }
      }
      if(type=="rstpm"){
        Stihat <- predict(model,
                          newdata=data.frame("eventtime"=rep(eventtime[cases[1]], length(comparable)+1),
                                             X[c(cases[i],comparable), ]),
                          type="surv")
      }
      conc[i, comparable] <- ifelse(Stihat[1] < Stihat[-1], 1, 0)
      conc[i, comparable] <- ifelse(Stihat[1] == Stihat[-1], .5, conc[i, comparable])
    }
  }
  list("Ctd"=sum(conc)/sum(comp),
       "propComparable"=sum(comp) / (length(eventtime)^2 - length(eventtime)))
}

# Time-fixed c-statistic
ctd.fixed <- function(eventtime, status, times, Shat){
  cases <- which(status == 1)
  comp <- matrix(0, length(eventtime), length(eventtime))
  conc <- array(0, dim=c(length(eventtime), length(eventtime), length(times)))
  for(i in 1:length(cases)){
    comparable <- which(eventtime > eventtime[cases[i]] | (eventtime == eventtime[cases[i]] & status %in% 0))
    if(length(comparable) > 0){
      comp[i, comparable] <- 1
      for(j in 1:length(times)){
        conc[i, comparable, j] <- ifelse(Shat[cases[i],j] < Shat[comparable,j], 1, 0)
        conc[i, comparable, j] <- ifelse(Shat[cases[i],j] == Shat[comparable,j], .5, conc[i, comparable, j])
      }
    }
  }
  list("Ctd"=sapply(1:length(times), function(x) sum(conc[,,x])/sum(comp)),
       "propComparable"=sum(comp) / (length(eventtime)^2 - length(eventtime)))
}

# clean-up to avoid running out of memory and keep track of any warnings
cleanup <- function(ws.names, print=FALSE){
  rm(list=ls(envir = globalenv())[!ls(envir = globalenv()) %in% ws.names], envir = globalenv())
  w <- warnings()
  if(!is.null(w)){
    assign("last.warning", NULL, envir = baseenv())
  }
  gc()
  if(print) return(w)
}
