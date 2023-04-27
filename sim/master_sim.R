# Master script to run the simulation study

# (NOT RUN) parts are commented out to avoid accidental runs that take a lot of time
# NB. Replicating the entire simulation study on a single machine would not be a reasonable endeavor.

# # Main simulations (NOT RUN)
# nsim <- 1000
# for(iter in 1:nsim){
#   set.seed(iter)
#   
#   sample.size.settings <- 1:4 # n=100, 250, 500, 1000
#   highercens <- FALSE
#   higherp <- FALSE 
#   
#   source("sim/sims.R")
#   
#   save(results, file= paste0("intermediates/", "survsim", iter, ".RData"))
# }


# # Simulation with a higher censoring rate (NOT RUN)
# nsim <- 500
# for(iter in 1:nsim){
#   set.seed(iter)
#   
#   sample.size.settings <- 2:4  # n=100, 250, 500, 1000
#   highercens <- TRUE
#   higherp <- FALSE 
#   
#   source("sim/sims.R")
#   
#   save(results, file= paste0("intermediates/", "survsim_highcens", iter, ".RData"))
# }

# # Higher number of covariates (NOT RUN)
# nsim <- 500
# for(iter in 1:nsim){
#   set.seed(iter)
#   
#   sample.size.settings <- 2:4 # n=100, 250, 500, 1000
#   highercens <- FALSE
#   higherp <- TRUE 
#   
#   source("sim/sims.R")
#   
#   save(results, file= paste0("intermediates/", "survsim_highp", iter, ".RData"))
# }

# iter 10 check MBP versus HPC
load("intermediates/otherMachine/survsim10.RData") # other machine replicate N=100
v1 <- results
load("intermediates/survsim10.RData") # HPC version
v2 <- results
v1$rmse[ ,, 1] 
v2$rmse[ ,, 1] 
all.equal(v1$rmse[1,,1],v2$rmse[1,,1])
all.equal(v1$cstats[[1]][[1]],v2$cstats[[1]][[1]])
