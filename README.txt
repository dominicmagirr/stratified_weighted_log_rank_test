## R version 4.1.0 (2021-05-18)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux
## 
## Matrix products: default
## BLAS/LAPACK: /CHBS/apps/EB/software/imkl/2019.1.144-gompi-2019a/compilers_and_libraries_2019.1.144/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] compiler_4.1.0  magrittr_2.0.3  fastmap_1.1.0   cli_3.4.1      
##  [5] tools_4.1.0     htmltools_0.5.3 rstudioapi_0.14 yaml_2.3.6     
##  [9] stringi_1.7.8   rmarkdown_2.17  knitr_1.40      stringr_1.4.1  
## [13] xfun_0.34       digest_0.6.30   rlang_1.0.6     evaluate_0.17


-------------------------------------------------------------
-- STEPS TO REPRODUCE THE OUTPUTS PRESENTED IN THE ARTICLE --
-------------------------------------------------------------

-- STEP 1: DOWNLOAD OAK AND POPLAR TRIAL DATA:

Go to the paper: https://doi.org/10.1038/s41591-018-0134-3, section 
"Supplementary Information", and then click on "Supplementary Table 8" to
(automatically) download the file "41591_2018_134_MOESM3_ESM.xlsx". Place 
this file in the folder "Data".

-- STEP 2: PRODUCE FIGURES 1 AND 2 (AND TABLE 2):

To reproduce Figures 1 and 2, and Table 2 (OAK and POPLAR results), you first to complete
STEP 1. Then you need to run the file "OAK_POPLAR_trials.R".

-- STEP 3: PRODUCE FIGURE 3:

To reproduce Figure 3 (visualization of weights), run the file "plot_weights.R".

-- STEP 4: PRODUCE FIGURES 4 AND 5:

To reproduce the main simulation results (Figure 4 and Figure 5), first run the file "simulate_scenarios_A.R". Because running this program takes some time, we leave as .RData files the results from the simulated scenarios, one using s_star and the other one using t_star. Next run "simulate_scenarios_B.R". This takes the saved .RData files and produces the figures.

-- STEP 5: PRODUCE FIGURES 6, 7 AND 8

To reproduce Figures 6,7,8 (i.e., the simulation scenarios) run the file "plot_scenarios.R".
