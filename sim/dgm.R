## sims ##
library(regsurv)
library(CVXR)
library(simsurv)
library(survival)
library(MASS)
library(glmnet)
library(rstpm2)

## Final DGM #########################################
set.seed(123)

S0 <- function(lambdas, gammas, pmix, t){
  (pmix * exp(-lambdas[1] * t^(gammas[1])) + (1-pmix) * exp(-lambdas[2] * t^(gammas[2])))
}
H0 <- function(lambdas, gammas, pmix, t){
  -log(pmix * exp(-lambdas[1] * t^(gammas[1])) + (1-pmix) * exp(-lambdas[2] * t^(gammas[2])))
}
h0 <- function(lambdas, gammas, pmix, t){
  (lambdas[1] * gammas[1] * t^(gammas[1]-1) * pmix * exp(-lambdas[1] * t^(gammas[1])) +
     lambdas[2] * gammas[2] * t^(gammas[2]-1) * (1-pmix) * exp(-lambdas[2] * t^(gammas[2]))) /
    (pmix * exp(-lambdas[1] * t^(gammas[1])) + (1-pmix) * exp(-lambdas[2] * t^(gammas[2])))
}
truehazard <- function(eventtime, X, groundtruth){
  with(groundtruth, h0(lambdas, gammas, pmix, t=eventtime) *
         as.vector(exp(as.matrix(X[, names(betas)]) %*% betas +
                         as.matrix(X[, names(tde)] * tdefunction(eventtime)) %*% tde)))
}
truecumhazard <- function(eventtime, X, groundtruth){
  sapply(1:nrow(X), function(s){
    stats::integrate(truehazard, lower=0, upper=eventtime[s],
                     X=X[s, ], groundtruth=groundtruth)$value
  })
}

rm(list=ls()[!ls() %in% c("h0","H0","S0","truecumhazard", "truehazard")])

maxt <- 30
t <- seq(0, maxt, 0.01)
lambdas <- c(.25,.6)/12
gammas <- c(1.1,1.4)
pmix <- 0.4
betas = c("x1" = 0, "x2"= 0, "x3"=.5, "x4"=-.5, "x5"=.25, "x6"=-.25, "x7"=.125, "x8"=-.125, "x9"=0.0625, "x10"=-0.0625, "x11"=.5)
tdefunction <- function(t) .9^t
tde <- c("x1" = -1, "x2"=.75, "x3"=-.5)

# pdf(file="figures/dgm.pdf",
#     width = .8*8, height = .8*6.6)
setEPS()
postscript(file="figures/dgm.eps",
           width = .8*8, height = .8*6.6)
par(mfrow=c(2,3))
plot(S0(lambdas, gammas, pmix, t) ~ t, type="l", ylim=c(0,1),
     main=expression(paste("S"[0])), ylab="Survival probability", xlab="Time")
plot(H0(lambdas, gammas, pmix, t) ~ t, type="l",
     main=expression(paste("H"[0])), ylab="Cumulative hazard", xlab="Time")
plot(h0(lambdas, gammas, pmix, t) ~ t, type="l",
     main=expression(paste("h"[0])), ylab="Hazard", xlab="Time")

# par(mfrow=c(2,3))
plot(betas[1] + tde[1] * tdefunction(t) ~ t, type="l", ylim=c(-1,1),
     main=expression(paste(beta[1](t))), ylab="log hazard ratio", xlab="Time"); abline(h=0, lty=3)
plot(betas[2] + tde[2] * tdefunction(t) ~ t, type="l", ylim=c(-1,1),
     main=expression(paste(beta[2](t))), ylab="log hazard ratio", xlab="Time"); abline(h=0, lty=3)
plot(betas[3] + tde[3] * tdefunction(t) ~ t, type="l", ylim=c(-1,1),
     main=expression(paste(beta[3](t))), ylab="log hazard ratio", xlab="Time"); abline(h=0, lty=3)
dev.off()

# # on log time scale (not required for the manuscript)
# plot(betas[1] + tde[1] * tdefunction(t) ~ log(t), type="l", ylim=c(-1,1)); abline(h=0, lty=3)
# plot(betas[2] + tde[2] * tdefunction(t) ~ log(t), type="l", ylim=c(-1,1)); abline(h=0, lty=3)
# plot(betas[3] + tde[3] * tdefunction(t) ~ log(t), type="l", ylim=c(-1,1)); abline(h=0, lty=3)

p <- length(betas)
rho <- .25
R <- matrix(rho, p, p)
diag(R) <- 1
D <- rep(1, p)
Sigma <- D * R * D

# # NOT RUN (takes 27 minutes on a 2022 MBP with Apple M1 pro chip)
# npop <- 110000
# X <- mvrnorm(n=npop, mu=rep(0, p), Sigma = Sigma)
# colnames(X) <- paste0("x", 1:p)
# ptm <- proc.time()
# population <- simsurv(dist = "weibull", lambdas = lambdas, gammas = gammas,
#                    mixture = TRUE, pmix=pmix, tde=tde, tdefunction = tdefunction,
#                    x = data.frame(X[ ,grep("x", colnames(X))]),
#                    betas=betas, maxt = maxt)
# population <- cbind(population, X)
# (proc.time() - ptm) / 3600
# validation <- population[100001:110000, ] # 10.000 cases for validation purposes
# dev_population <- population[-c(100001:110000), ] # the remainder to sample development sets from.
# save(validation, file="data/validation.RData")
# save(dev_population, file="data/dev_population.RData")
# 
# load(file="data/validation.RData")
# load(file="data/dev_population.RData")
# 
# # add noise variables for higherp simulations
# set.seed(456)
# head(dev_population)
# 
# noise <- matrix(rnorm(nrow(dev_population) * 20), nrow(dev_population), 20)
# colnames(noise) <- paste0("n", 1:20)
# dev_population <- cbind(dev_population, noise)
# 
# noise <- matrix(rnorm(nrow(validation) * 20), nrow(validation), 20)
# colnames(noise) <- paste0("n", 1:20)
# validation <- cbind(validation, noise)
# 
# save(validation, file="data/validation.RData")
# save(dev_population, file="data/dev_population.RData")

load(file="data/validation.RData")
load(file="data/dev_population.RData")

# ## Checks on DGM (not required for the manuscript)
# # Split in several instances due to very large data set being constructed when using tt() (start stop splits for all event times)
# check.dgm <- list()
# index <- 1:1000
# for(a in 0:9){
#   check.dgm[[a+1]] <- coxph(Surv(eventtime, status) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+tt(x1)+tt(x2)+tt(x3),
#                             tt=function(x,t,...) x*tdefunction(t), data=validation[index+a*1000, ])
#   print(a+1)
# }
# apply(sapply(check.dgm, coef), 1, mean) # check
# boxplot(t(sapply(check.dgm, coef))) # check
# abline(h=c(betas, tde))
# rm(check.dgm, index, a)

# Store ground truth (true survival probabilities, hazards and cumulative hazard) for the validation data
# # (NOT RUN) takes 6 minutes
# ptm <- proc.time()
# groundtruth <- list()
# groundtruth$betas <- betas
# groundtruth$pmix <- pmix
# groundtruth$lambdas <- lambdas
# groundtruth$gammas <- gammas
# groundtruth$tde <- tde
# groundtruth$tdefunction <- tdefunction
# groundtruth$maxt <- maxt
# groundtruth$n <- npop
# groundtruth$p <- p
# groundtruth$covmu <- rep(0, p)
# groundtruth$covSigma <- Sigma
# rm(R, D, Sigma, X, p, betas, lambdas, gammas, pmix, t, rho, tde, tdefunction)
# 
# # hazard for the observed event times
# htrue <- truehazard(eventtime=validation$eventtime,
#                     X=as.matrix(validation[ ,c(grep("x", names(validation)))]),
#                     groundtruth=groundtruth)
# groundtruth$htrue <- htrue
# 
# # cumulative hazard for the observed event times
# Htrue <- truecumhazard(eventtime=validation$eventtime,
#                        X=validation[ ,grep("x", names(validation))],
#                        groundtruth=groundtruth)
# groundtruth$Htrue <- Htrue
# 
# # cumulative hazard the fixed time points of interest,
# # being 2.5, 5, 7.5, 10, 20, and 30t
# Htrue2.5 <- truecumhazard(eventtime=rep(2.5, length(validation$eventtime)),
#                           X=validation[ ,grep("x", names(validation))],
#                           groundtruth=groundtruth)
# groundtruth$Htrue2.5 <- Htrue2.5
# Htrue5 <- truecumhazard(eventtime=rep(5, length(validation$eventtime)),
#                         X=validation[ ,grep("x", names(validation))],
#                         groundtruth=groundtruth)
# groundtruth$Htrue5 <- Htrue5
# Htrue7.5 <- truecumhazard(eventtime=rep(7.5, length(validation$eventtime)),
#                           X=validation[ ,grep("x", names(validation))],
#                           groundtruth=groundtruth)
# groundtruth$Htrue7.5 <- Htrue7.5
# Htrue10 <- truecumhazard(eventtime=rep(10, length(validation$eventtime)),
#                          X=validation[ ,grep("x", names(validation))],
#                          groundtruth=groundtruth)
# groundtruth$Htrue10 <- Htrue10
# Htrue20 <- truecumhazard(eventtime=rep(20, length(validation$eventtime)),
#                          X=validation[ ,grep("x", names(validation))],
#                          groundtruth=groundtruth)
# groundtruth$Htrue20 <- Htrue20
# Htrue30 <- truecumhazard(eventtime=rep(30, length(validation$eventtime)),
#                          X=validation[ ,grep("x", names(validation))],
#                          groundtruth=groundtruth)
# groundtruth$Htrue30 <- Htrue30
# save(groundtruth, file="data/groundtruth.RData")
# (proc.time() - ptm) / 3600
load("data/groundtruth.RData")

# (not required for the manuscript)
# Visualize true hazard, cumulative hazards and survival probs at observed event times
par(mfrow=c(1,3))
s <- sample(1:nrow(validation), 1000)
plot(log(groundtruth$htrue)[s] ~ validation$eventtime[s], main="true log hazard")
mtext("log(h0) in red")
lines(log(with(groundtruth, h0(lambdas, gammas, pmix, t=seq(0,maxt,.1)))) ~
        seq(0,maxt,.1), col="red")

plot(log(groundtruth$Htrue)[s] ~ validation$eventtime[s], main="true log hazard")
mtext("log(H0) in red")
lines(log(with(groundtruth, H0(lambdas, gammas, pmix, t=seq(0,maxt,.1)))) ~
        seq(0,maxt,.1), col="red")

Strue <- exp(-groundtruth$Htrue)
plot(Strue[s] ~ validation$eventtime[s], main="true log hazard")
mtext("S0 in red")
lines(with(groundtruth, S0(lambdas, gammas, pmix, t=seq(0,maxt,.1))) ~
        seq(0,maxt,.1), col="red")