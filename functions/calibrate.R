#' Calibration according to Austin 2020 (https://doi.org/10.1002/sim.8570)
#'
#' @param Shat Predicted survival probabilities for the particular time point
#' @param time The time point corresponding to the prediction time for Shat
#' @param validation # validation data.frame that should have colums 'eventtime' and 'status' (0 for censored, 1 for events)
#' @param predict.grid 
#' @param plot FALSE for silent output; TRUE for the calibration plot
#' @param store.predict.grid TRUE if the predict.grid input should be stored in the output object
#'
#' @return a list containing all relevant output
#' \item{cox.cal.stats}{contains the ICI, E50 and E90 based on CoxPH calibration and a test for non-
#' proportional hazards (indicating the need for hazard regression calibration)}
#' \item{cox.cal.obs}{stores the modeled ('observed') survival probabilities along predict.grid,
#' which can be used to construct a calibration plot}
#' \item{hare.cal.stats}{stores the ICI, E50 and E90 based on hazard regression calibration}
#' \item{hare.cal.obs}{stores the modeled ('observed') survival probabilities along predict.grid,
#' which can be used to construct a calibration plot}
#' \item{predict.grid}{return the predict.grid input if predict.grid=TRUE}
#' \itme{plots}{if plot=TRUE, returns a plot for both CoxPH and hazard regression calibration with
#' the calibration curve in solid red, the density of predicted probabilities in dotted black,
#' and the reference line for perfect calibration in dashed black}
#' 
#' @examples
#' library(survival)
#' names(veteran)[which(names(veteran) == "time")] <- "eventtime"
#' mod <- coxph(Surv(eventtime, status) ~ trt + celltype + karno +
#'                 diagtime + age + prior, data=veteran)
#' newdata <- veteran
#' newdata$eventtime <- 60 # predict for eventtime 60 days                 
#' Shat <- exp(-predict(mod, newdata=newdata, type="expected"))
#' 
#' # the example just uses the same veteran data, but validation data should in practice be
#' # out-of-sample data (bootstrap, CV, external data)
#' calibrate(Shat=Shat, time=60, validation=veteran, plot=TRUE)
#' 
#' # data split example
#' library(survival)
#' # The rotterdam data set includes 2982 primary breast cancers patients whose
#' # records were included in the Rotterdam tumor bank.
#' rotterdam$status  <- pmax(rotterdam$recur, rotterdam$death)
#' rotterdam$eventtime <- with(rotterdam, ifelse(recur==1, rtime, dtime))
#' ntrain <- 1750
#' set.seed(1)
#' train.index <- sample(1:nrow(rotterdam), ntrain)
#' train <- rotterdam[train.index, ]
#' mod <- coxph(Surv(eventtime, status) ~ pspline(age, df=4) + meno + size +
#'                er + pspline(nodes, df=4), data=train)
#' test <- rotterdam[-train.index, ]
#' newdata <- test
#' newdata$eventtime <- 5 * 365.25 # predict survival probabilites at 5 years
#' Shat <- exp(-predict(mod, newdata=newdata, type="expected"))
#' calibrate(Shat=Shat, time=1826, validation=test, plot=TRUE)
#' 
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
           xlab = "Predicted survival probability",
           ylab = "Observed survival probability",
           main = "Calibration (based on CoxPH)")
      abline(0,1,lty=2)
      par(new=T)
      plot(density(Fhat),axes=F,xlab=NA,ylab=NA,main="",xlim=c(0,1), lty=3)
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
         xlab = "Predicted survival probability",
         ylab = "Observed survival probability",
         main = "Calibration (based on hazard regression)")
    abline(0,1,lty=2)
    par(new=T)
    plot(density(Fhat),axes=F,xlab=NA,ylab=NA,main="",xlim=c(0,1), lty=3)
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

