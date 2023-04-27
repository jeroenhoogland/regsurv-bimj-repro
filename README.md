Reproducible research files
================

This README documents the procedures to replicate all analyses in the
manuscript entitled “Regularized parametric survival modeling to improve
risk prediction models”.

Authors: J Hoogland, TPA Debray, MJ Crowther, RD Riley, J IntHout, JB
Reitsma, AH Zwinderman Code was written by J Hoogland In case of
questions or comments please contact <j.hoogland@amsterdamumc.nl>

# Software versions

The code was written/evaluated in R with the following software
versions:

    #> R version 4.2.0 (2022-04-22)
    #> Platform: x86_64-apple-darwin17.0 (64-bit)
    #> Running under: macOS Big Sur/Monterey 10.16
    #> 
    #> Matrix products: default
    #> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    #> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    #> 
    #> locale:
    #> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    #> 
    #> attached base packages:
    #> [1] splines   stats     graphics  grDevices utils     datasets  methods  
    #> [8] base     
    #> 
    #> other attached packages:
    #>  [1] stringr_1.5.0      Hmisc_4.7-2        ggplot2_3.4.0      Formula_1.2-4     
    #>  [5] lattice_0.20-45    polspline_1.1.22   rstpm2_1.5.8       glmnet_4.1-6      
    #>  [9] Matrix_1.5-3       MASS_7.3-58.1      survival_3.4-0     simsurv_1.0.0     
    #> [13] CVXR_1.0-11        regsurv_0.0.0.9000 dplyr_1.0.10      
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] Rcpp_1.0.10         bdsmatrix_1.3-6     mvtnorm_1.1-3      
    #>  [4] deldir_1.0-6        png_0.1-8           assertthat_0.2.1   
    #>  [7] digest_0.6.31       foreach_1.5.2       utf8_1.2.2         
    #> [10] gmp_0.7-1           R6_2.5.1            backports_1.4.1    
    #> [13] stats4_4.2.0        evaluate_0.19       pillar_1.8.1       
    #> [16] rlang_1.0.6         data.table_1.14.6   rstudioapi_0.14    
    #> [19] rpart_4.1.19        checkmate_2.1.0     bbmle_1.0.25       
    #> [22] rmarkdown_2.19      foreign_0.8-84      htmlwidgets_1.6.0  
    #> [25] bit_4.0.5           munsell_0.5.0       compiler_4.2.0     
    #> [28] numDeriv_2016.8-1.1 xfun_0.36           pkgconfig_2.0.3    
    #> [31] base64enc_0.1-3     shape_1.4.6         mgcv_1.8-41        
    #> [34] htmltools_0.5.4     nnet_7.3-18         tidyselect_1.2.0   
    #> [37] htmlTable_2.4.1     gridExtra_2.3       tibble_3.1.8       
    #> [40] codetools_0.2-18    fansi_1.0.3         withr_2.5.0        
    #> [43] grid_4.2.0          nlme_3.1-161        gtable_0.3.1       
    #> [46] lifecycle_1.0.3     DBI_1.1.3           magrittr_2.0.3     
    #> [49] scales_1.2.1        cli_3.5.0           stringi_1.7.8      
    #> [52] Rmpfr_0.9-1         latticeExtra_0.6-30 generics_0.1.3     
    #> [55] vctrs_0.5.1         deSolve_1.34        RColorBrewer_1.1-3 
    #> [58] iterators_1.0.14    tools_4.2.0         interp_1.1-3       
    #> [61] bit64_4.0.5         glue_1.6.2          jpeg_0.1-10        
    #> [64] fastmap_1.1.0       yaml_2.3.6          colorspace_2.0-3   
    #> [67] cluster_2.1.4       knitr_1.41

As an exception to the above, simulations were run in parallel on a
Linux server with software versions:

    #> R version 3.6.3 (2020-02-29)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: CentOS Linux 7 (Core)
    #> 
    #> Matrix products: default
    #> BLAS:   /hpc/local/CentOS7/julius_te/R-3.6.3/lib64/R/lib/libRblas.so
    #> LAPACK: /hpc/local/CentOS7/julius_te/R-3.6.3/lib64/R/lib/libRlapack.so
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> attached base packages:
    #> [1] splines   stats     graphics  grDevices utils     datasets  methods  
    #> [8] base     
    #> 
    #> other attached packages:
    #>  [1] Hmisc_4.7-0        ggplot2_3.3.5      Formula_1.2-4      lattice_0.20-38   
    #>  [5] polspline_1.1.19   rstpm2_1.5.2       glmnet_4.1-2       Matrix_1.2-18     
    #>  [9] MASS_7.3-51.5      survival_3.3-1     simsurv_1.0.0      CVXR_1.0-9        
    #> [13] regsurv_0.0.0.9000 dplyr_1.0.7       
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] Rcpp_1.0.7          bdsmatrix_1.3-4     mvtnorm_1.1-2      
    #>  [4] png_0.1-7           digest_0.6.27       foreach_1.5.1      
    #>  [7] utf8_1.2.1          gmp_0.6-2           R6_2.5.1           
    #> [10] backports_1.3.0     stats4_3.6.3        pillar_1.6.1       
    #> [13] rlang_0.4.11        data.table_1.14.2   rstudioapi_0.13    
    #> [16] rpart_4.1-15        bbmle_1.0.23.1      checkmate_2.0.0    
    #> [19] stringr_1.4.0       foreign_0.8-75      htmlwidgets_1.5.4  
    #> [22] bit_4.0.4           munsell_0.5.0       xfun_0.24          
    #> [25] compiler_3.6.3      numDeriv_2016.8-1.1 pkgconfig_2.0.3    
    #> [28] base64enc_0.1-3     shape_1.4.6         mgcv_1.8-31        
    #> [31] htmltools_0.5.2     nnet_7.3-12         tidyselect_1.1.1   
    #> [34] tibble_3.1.2        gridExtra_2.3       htmlTable_2.3.0    
    #> [37] codetools_0.2-16    fansi_0.5.0         crayon_1.4.1       
    #> [40] withr_2.4.2         grid_3.6.3          orthopolynom_1.0-5 
    #> [43] nlme_3.1-144        gtable_0.3.0        lifecycle_1.0.0    
    #> [46] magrittr_2.0.1      scales_1.1.1        stringi_1.6.2      
    #> [49] Rmpfr_0.8-7         latticeExtra_0.6-29 ellipsis_0.3.2     
    #> [52] generics_0.1.1      vctrs_0.3.8         deSolve_1.28       
    #> [55] RColorBrewer_1.1-2  tools_3.6.3         iterators_1.0.13   
    #> [58] bit64_4.0.5         glue_1.4.2          purrr_0.3.4        
    #> [61] jpeg_0.1-9          fastmap_1.1.0       colorspace_2.0-2   
    #> [64] cluster_2.1.0       ECOSolveR_0.5.4     knitr_1.33

# File structure

- `README.html` the general README information you’re now reading.
- Folder `Results` contains
  - `Results_bimj.Rmd` contains the Markdown script to generate
    `Results_bimj.html`, which generates all results for the manuscript.
- Folder `./sim/` contains all R scripts pertaining to the simulation
  study
  - `dgm.R` contain the R script to generate all simulated data (files
    `groundtruth.RData`, `dev_population.RData`, and `validation.RData`)
    and the supplementary in the data generating mechanims section
    (`dgm.eps`).
  - `master_sim.R` Master file that runs all simulations and stores the
    results. Note that each run has its own random seed, such that
    individual runs can be replicated. NB: All sample size settings are
    generated from the same run.
  - `groundtruth.RData` contains the ground truth for the simulated
    data.
  - `dev_population.RData`contains the simulated development population
    (N=100,000)
  - `validation.RData` contains the simulated validation data
  - `helperFunctions.R` contains helper function that are not packaged
    (only useful to avoid code replication in the current context)
  - folder `./figures/` contain all manuscript figures in the required
    formats
    <!-- * Folder `./HPC_bimj/` contains scripts for the parallel runs on a linux server (not required) -->
    <!-- * `*.sh` shell scripts that contains the code to run the simulations on a CentOS High Performance Grid (HPC) using Slurm. -->
    <!-- * `*.R` contains the R scripts for a single run for each group of settings (main, high censoring, high p, applied example). -->
- Folder `./ae/` contains all R scripts pertaining to the applied
  example
  - `ae_apparent.R` contains the R script for the main analyses
  - `ae_bimj_boot.R` contains the R script for a single bootstrap run
  - `helperFunctions_ae.R` contains helper function that are not
    packaged (only useful to avoid code replication in the current
    context).
  - `ae_apparent_bimj.RData` contains the pre-computed results of
    running `ae_apparent.RData`.
- Folder `./functions` contains
  - The documented function `calibrate()` stored in `calibrate.R`, which
    can be loaded using `source("functions/calibrate.R")`
  - (Note that all analysis functions are directly loaded from R
    packages)
- Folder `./intermediates/` contains all intermediate files (output
  files) from the high performance grid for both the simulations
  (`survsim*.RData`, `survsim_highp*.RData`, `survsim_highcens*.RData`)
  and the applied example (`ae*.RData`), where `*` is replaced by the
  seed number for the random number generator.

# Reproduce all figures and numerical results

The required computation time for the simulations is long (on average
nearly 5 hours for the main simulation runs on the CentOS system and
longer for the additional higher dimensional and increased censoring
settings).

Therefore, intermediate files are available for all simulation runs and
the bootstrap runs for the applied example.

## The Results_bimj markdown and html

`Results_bimj.Rmd` contains all code to constructs all figures and
numerical results from the intermediate files, and the corresponding
`html` file `Results_bimj.html` provides the rendered output. Note that
numerical results are provided in `Results_bimj.html` (such as Table 1),
and that figures are not rendered but stored in the correct format in
folder `./figures`. Each of the figures has the exact same file name as
in the Latex files, such that they can be inserted in the right place
without any modification.

## Replicating the simulation data

R script `dgm.R` in folder `./sim/` generates all required simulated
data and stored it in the right place. Also, it generates `dgm.eps`,
which is the supplementary figure on the data generating mechanism.

Resulting files are stored in folder `./data/` and are named
`dev_population.RData` (the part of the simulated population from which
development data sets / training data are sampled), `validation.RData`
(the part of the simulated population reserved for validation purposes),
and `groundtruth.RData` (which contains all necessary information on the
true model and true survival probabilities).

## Replicating the intermediate files

First of all, R package **regsurv** needs to be installed from GitHub

    library(devtools)
    install_github("jeroenhoogland/regsurv")

### Simulations

In folder `./sim/`, R script `master_sim.R` can run all simulations (for
all simulation settings). The specifications for the different
simulation settings can be set within the master file and have intuitive
names: `sample.size.settings` takes number 1:4 for each of the 4 sample
size settings in the manuscript (N=100, 250, 500, 1000), `highercens` is
TRUE for the high censoring settings, `higherp` is true for the setting
with a larger set of covariates, and both are FALSE for the default
simulation settings.

An impression of the script to run the main simulation settings is

    # Master script to run the simulation study
    nsim <- 1000

    # Main simulations (NOT RUN)
    for(iter in 1:nsim){
      set.seed(iter)
      
      sample.size.settings <- 1:4 # n=100, 250, 500, 1000
      highercens <- FALSE
      higherp <- FALSE 
      
      source("sim/sims.R")
      
      save(results, file= paste0("intermediates/", "survsim", iter, ".RData"))
    }

R script `master_sim.R` also includes a few lines of code that show spot
replication of the 10th simulation run for sample size N=100 on a
different machine.

As a generic example, to replicate let’s say run 98 for sample size
setting 1 of the main simulations, the following code can be used:

    iter <- 98
    set.seed(iter)
    sample.size.settings <- 1
    highercens <- FALSE
    higherp <- FALSE
    source("sim/sims.R")
    spot_check <- results

On a recent MBP (2021, Apple M1 Pro, 16GB working memory), this takes
approximately 48 minutes.

This generates a results object that can be checked against
`./intermediates/survsim98.RData`.

For instance, to check the overall root mean squared prediction errors
for each of the methods, compare

    load("intermediates/survsim98.RData") # loads the results object for sim 98
    all.equal(spot_check$rmse[1,,], results$rmse[1,,], tolerance=10^-5) # should return TRUE
    all.equal(sapply(spot_check$cstats[[1]], function(x) x[1, ]), 
              sapply(results$cstats[[1]], function(x) x[1, ]), tolerance=10^-5) # should return TRUE

### Applied example

The intermediate files for the applied example can also be reproduced.
The apparent (whole sample) results can be replicated using R script
`ae_apparent.R` in the `./ae/` folder. The bootstrap runs can be
replicated using `ae_bimj_boot.R` in the same folder.

NB boot needs b and seed

For instance:

    b <- 9 # for replication of bootstrap 9
    set.seed(b)
    source("ae/ae_bimj_boot.R")
    spot_check <- boot

Checks against

    load("intermediates/ae9.RData") 
    unlist(spot_check$cstats.oob) # OOB Ctd estimates
    all.equal(unlist(spot_check$cstats.oob), unlist(boot$cstats.oob), tol=10^-5)

END
