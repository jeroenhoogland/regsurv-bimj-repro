# generic plot function for time, rmspe, time-dependent rmspe, and discrimination results
plot.array.results <- function(x, main, ylab, ylim=NULL, plot.min=FALSE, plot.max=FALSE){
  boxplot(t(x), las=1, ylab=ylab, xlab="", main=main, ylim=ylim, xaxt="n")
  axis(1, las=2, at = 1:10, labels = methods)
  mm <- apply(x, 1, function(x) median(x, na.rm=TRUE))
  if(plot.min) abline(h=mm[mm == min(mm)], lty=3)
  if(plot.max) abline(h=mm[mm == max(mm)], lty=3)
}

# extract processing times
times <- array(sapply(sr$time.per.model.sec, function(x) x), dim=c(nn,nmethods, nsim))

if(filetype == ".eps"){
  setEPS()
  postscript(file=paste0(dir, "procTimes", postfix, filetype),
             width = 10, height = 10)
} else if (filetype==".pdf"){
  pdf(file=paste0(dir, "procTimes", postfix, filetype),
      width = 10, height = 10)
}
par(mfrow=c(2,2))
par(font.main = 1)
# par(mar = c(bottom, left, top, right)) # set margins; default is c(5.1, 4.1, 4.1, 2.1)
par(mar = c(7.1, 4.1, 4.1, 2.1))
for(i in sample.size.settings){
  plot.array.results(log(times[i,,]), 
                     main=c("N=100", "N=250","N=500", "N=1000")[i],
                     ylab="Log(time in seconds)", ylim=c(-5.4, 9))
}
dev.off()

# rmspe
rmspe <- array(sapply(sr$rmse, function(x) x[,,"overall"]), dim=c(nn,nmethods, nsim))
ylim <- list(c(0.12, 0.32), c(0.12, 0.24), c(0.12, 0.24), c(0.12, 0.24))
if(filetype == ".eps"){
  setEPS()
  postscript(file=paste0(dir, "rMSPE", postfix, filetype),
             width = 10, height = 10)
} else if (filetype==".pdf"){
  pdf(file=paste0(dir, "rMSPE", postfix, filetype),
      width = 10, height = 10)
}
par(mfrow=c(2,2))
par(font.main = 1)
par(mar = c(7.1, 4.1, 4.1, 2.1))
for(i in sample.size.settings){
  plot.array.results(rmspe[i,,], 
                     main=c("N=100", "N=250","N=500", "N=1000")[i],
                     ylab="rMSPE", ylim=ylim[[i]],
                     plot.min = TRUE)
}
dev.off()

# rmspe for specific time points
ntimes <- 6
rmspet <- array(sapply(sr$rmse, function(x) x[,,-1]), dim=c(nn,nmethods, ntimes, nsim))
dimnames(rmspet) <- list(c("N=100", "N=250","N=500", "N=1000"),
                         colnames(sr$rmse$sim1),
                         c("t=2.5", "t=5", "t=7.5", "t=10", "t=20", "t=30"),
                         names(sr$rmse))
TV <- c(1,2,4,6,8,10)
PH <- (1:10)[-TV]
col <- RColorBrewer::brewer.pal(10, "Paired")
if(filetype == ".eps"){
  setEPS()
  postscript(file=paste0(dir, "rMSPEt", postfix, filetype),
             width = 10, height = 10)
} else if (filetype==".pdf"){
  pdf(file=paste0(dir, "rMSPEt", postfix, filetype),
      width = 10, height = 10)
}
par(font.main = 1)
par(mfrow=c(2,2))
if(1 %in% sample.size.settings){
  ylim1 <- c(min(t(sapply(sample.size.settings, 
                          function(x) apply(rmspet[x,,,], c(1,2), function(x) mean(x, na.rm=TRUE))))), 
             max(apply(rmspet[1,,,], c(1,2), function(x) mean(x, na.rm=TRUE))))
  ylim234 <- (range(t(sapply(sample.size.settings[sample.size.settings %in% 2:4], function(x){
    apply(rmspet[x,,,], c(1,2), function(x) mean(x, na.rm=TRUE))}))) + ylim1) / 2
} else {
  ylim234 <- range(t(sapply(sample.size.settings[sample.size.settings %in% 2:4], function(x){
    apply(rmspet[x,,,], c(1,2), function(x) mean(x, na.rm=TRUE))})))
}
for(i in sample.size.settings){
  matplot(c(2.5,5,7.5,10,20,30), 
          t(apply(rmspet[i,TV,,], c(1,2), function(x) mean(x, na.rm=TRUE))), 
          type="l", lty=1, lwd=2, col=col[TV],
          main=c("N=100", "N=250","N=500", "N=1000")[i], ylab="rMSPE", xlab="Time (in months)",
          ylim=if(i == 1) ylim1 else ylim234)
  matplot(c(2.5,5,7.5,10,20,30), 
          t(apply(rmspet[i,PH,,], c(1,2), function(x) mean(x, na.rm=TRUE))), 
          type="l", lty=2, lwd=2, col=col[PH],
          add=TRUE)
  legend("bottomright",
         legend=methods,
         col=col, lty= ifelse(1:10 %in% TV,1,2), lwd=2, cex=1.1, ncol=2)
}
dev.off()

## Discrimination (supplementary)
filetype <- ".pdf" # (multi-page not supported for .eps)
# extract c-statistic results (overall and for fixed time-points)
cstat <- lapply(sr$cstats, function(x)
  array(sapply(x, function(xx) 
      sapply(xx, function(xxx) 
        if(is.null(xxx)) rep(NA, 7) else xxx["C Index", ])),
    dim=c(ntimes+1, nmethods, nn)))
cstat <- array(unlist(cstat), dim=c(ntimes+1,  nmethods, nn, nsim))
if(filetype == ".eps"){
  setEPS()
  postscript(file=paste0(dir, "cstat", postfix, filetype),
             width = 10, height = 10)
} else if (filetype==".pdf"){
  pdf(file=paste0(dir, "cstat", postfix, filetype),
      width = 10, height = 10)
}
tt <- c("(overall)", "(t=2.5)", "(t=5)", "(t=7.5)", "(t=10)", "(t=20)", "(t=30)")
par(mfrow=c(2,2))
par(font.main = 1)
par(mar = c(7.1, 4.1, 4.1, 2.1))
for(j in 1:7){
  for(i in sample.size.settings){
    plot.array.results(cstat[j,,i,], 
                       main=paste(c("N=100", "N=250","N=500", "N=1000")[i], tt[j]),
                       ylab="C-statistic", 
                       plot.max = TRUE)
    if(length(sample.size.settings)==3 & i==4) plot.new() 
  }
}
dev.off()

## Calibration (supplement)
if(filetype == ".eps"){
  setEPS()
  postscript(file=paste0(dir, "cal", postfix, filetype),
             width = 7, height = 10)
} else if (filetype==".pdf"){
  pdf(file=paste0(dir, "cal", postfix, filetype),
      width = 7, height = 10)
}
par(mfrow=c(3,2))
par(font.main = 1)
predict.grid <- 1-sr$predict.grid$sim1
plotCI90 <- FALSE
for(m in 1:10){
  for(t in 1:6){
    for(n in sample.size.settings){
      cal <- 1-sapply(sr$cal, function(x) 
        if(is.null(x[[n]][[m]][[t]]$hare.cal.obs)) rep(NA, nrow(predict.grid)) 
          else x[[n]][[m]][[t]]$hare.cal.obs)
      if(n==sample.size.settings[1]){
        plot(apply(cal, 1, function(x) mean(x, na.rm=TRUE)) ~ predict.grid[,t], 
                type="l", xlim=c(0,1), ylim=c(0,1), lwd=2,
                main=paste(methods[m],
                        paste0("t=",c(2.5,5,7.5,10,20,30))[t], sep = ", "),
             ylab="Observed", xlab="Predicted", col=n)
        legend("bottomright",
               legend=c("N=100", "N=250", "N=500", "N=1000")[sample.size.settings],
               col=sample.size.settings, lty=1, lwd=2, cex=.75, ncol=1)
        if(plotCI90){
          lines(apply(cal, 1, function(x) quantile(x, probs=.1, na.rm=TRUE)) ~ predict.grid[,t], lty=2, lwd=1)
          lines(apply(cal, 1, function(x) quantile(x, probs=.9, na.rm=TRUE)) ~ predict.grid[,t], lty=2, lwd=1)
        }
        abline(0,1,lty=1,col="red")
      } else {
        lines(apply(cal, 1, function(x) mean(x, na.rm=TRUE)) ~ predict.grid[,t], col=n, lwd=2)
        if(plotCI90){
          lines(apply(cal, 1, function(x) quantile(x, probs=.1, na.rm=TRUE)) ~ predict.grid[,t], lty=2, col=n, lwd=1)
          lines(apply(cal, 1, function(x) quantile(x, probs=.9, na.rm=TRUE)) ~ predict.grid[,t], lty=2, col=n, lwd=1)
        }
      }
    }
  }
}
dev.off()

