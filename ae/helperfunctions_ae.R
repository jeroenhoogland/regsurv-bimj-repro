# predict.cox.td.ae is model specific
predict.cox.td <- function(mod, newdata, ndtimes=NULL, type="expected", logtime=FALSE, iknots = NULL){
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
  
  pred.matrix <- survSplit(Surv(eventtime, status==1) ~., newdata, cut=Htimes, id="id")
  names(pred.matrix)[names(pred.matrix) == "event"] <- "status"
  if(logtime){
    pred.matrix$time.spline <- log(pred.matrix$eventtime)
  } else {
    pred.matrix$time.spline <- pred.matrix$eventtime
  }
  
  mmtt <- model.matrix( ~  id + tstart + eventtime + status + trt + celltype + karno +
                          diagtime + age + prior + (trt + celltype + karno +
                                                      diagtime + age + prior):time.spline,
                        data=pred.matrix)
  mmtt <- mmtt[ ,-1]
  vars <- c("trt", "celltypesmallcell", "celltypeadeno", "celltypelarge", "karno", "diagtime", "age", "prior")
  mmtt <- data.frame(mmtt)
  names(mmtt) <- c("id", "tstart", "eventtime", "status",vars,paste0(vars, "t1"))
  
  pred.matrix <- mmtt
  
  pred.matrix$H0hat <- predict(mod, newdata=pred.matrix, type="expected")
  H0hat.2 <- pred.matrix %>% group_by(id) %>%
    summarise(cumhaz=sum(H0hat))
  H0hat.2 <- as.numeric(H0hat.2$cumhaz)
  H0hat.2[is.na(H0hat.2)] <- 0
  H0hat.2
}
