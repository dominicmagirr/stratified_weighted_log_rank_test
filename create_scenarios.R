library(tidyverse)
library(survival)
library(ggpubr)
devtools::load_all("~/NPH/swlrt")

source("swlrt_fast.R")
source("wlrt_fast.R")

#######################################
surv_2 <- function(t, 
                   lambda_1, 
                   lambda_2, 
                   t_star){
  
  p_1 <- exp(-lambda_1 * t)
  p_2 <- exp(-lambda_1 * t_star) * exp(-lambda_2 * pmax(0, (t - t_star)))
  
  (t < t_star) * p_1 + (t >= t_star) * p_2
  
}

plot_strat <- function(t_end = 24, 
                       t_star_1 = 9,
                       t_star_0 = 9,
                       strat_1_rate_c_1 = log(2) / 4,
                       strat_1_rate_c_2 = log(2) / 4,
                       strat_1_rate_e_1 = log(2) / 4 * 0.6,
                       strat_1_rate_e_2 = log(2) / 4 * 0.6,
                       strat_0_rate_c_1 = log(2) / 15,
                       strat_0_rate_c_2 = log(2) / 15,
                       strat_0_rate_e_1 = log(2) / 15 * 1.1,
                       strat_0_rate_e_2 = log(2) / 15 * 1,
                       prop_strat_1 = 0.1,
                       label = "Non-prog.",
                       effect = "Null effect",
                       scenario = "Scenario 1"){

  t_seq <- seq(0, t_end, length.out = 100)
  
  s_c_1 <- surv_2(t_seq, lambda_1 = strat_1_rate_c_1, lambda_2 = strat_1_rate_c_2, t_star = t_star_1)
  s_e_1 <- surv_2(t_seq, lambda_1 = strat_1_rate_e_1, lambda_2 = strat_1_rate_e_2, t_star = t_star_1)
  s_c_0 <- surv_2(t_seq, lambda_1 = strat_0_rate_c_1, lambda_2 = strat_0_rate_c_2, t_star = t_star_0)
  s_e_0 <- surv_2(t_seq, lambda_1 = strat_0_rate_e_1, lambda_2 = strat_0_rate_e_2, t_star = t_star_0)
  

  
  first_stratum_df = data.frame(t_seq = rep(t_seq,2), 
                                surv = c(s_c_1,s_e_1),
                                arm = c(rep("Control",length(s_c_1)),rep("Experimental",length(s_e_1))),
                                stratum = "First stratum")
  
  second_stratum_df = data.frame(t_seq = rep(t_seq,2), 
                                surv = c(s_c_0,s_e_0),
                                arm = c(rep("Control",length(s_c_0)),rep("Experimental",length(s_e_0))),
                                stratum = "Second stratum")
  
  overall_df = data.frame(t_seq = rep(t_seq,2),
                          surv = c(prop_strat_1 * s_c_1 + (1 - prop_strat_1) * s_c_0, prop_strat_1 * s_e_1 + (1 - prop_strat_1) * s_e_0),
                          arm = c(rep("Control",length(s_c_1)),rep("Experimental",length(s_e_1))),
                          stratum = "Overall")
  
  plots_df = rbind(first_stratum_df, second_stratum_df, overall_df) %>%
    mutate(prognosis = label,
           effect = effect,
           scenario = scenario)
  

  
}



##########################################
## 1. Homogeneous
##########################################

###### A. Non-prognostic #################

plot_homogeneous_non_prognostic = plot_strat(t_end = 24, 
                                             t_star_1 = 6,
                                             t_star_0 = 6,
                                             strat_1_rate_c_1 = log(2) / 8,
                                             strat_1_rate_c_2 = log(2) / 8,
                                             strat_1_rate_e_1 = log(2) / 8 * (2/3),
                                             strat_1_rate_e_2 = log(2) / 8 * (2/3),
                                             strat_0_rate_c_1 = log(2) / 8,
                                             strat_0_rate_c_2 = log(2) / 8,
                                             strat_0_rate_e_1 = log(2) / 8 * (2/3),
                                             strat_0_rate_e_2 = log(2) / 8 * (2/3),
                                             prop_strat_1 = 0.5,
                                             label = "Non-prog.",
                                             effect = "PH",
                                             scenario = "Scenario 1")


###### B. moderate prognostic ############

plot_homogeneous_moderate_prognostic = plot_strat(t_end = 24, 
                                                  t_star_1 = 6,
                                                  t_star_0 = 6,
                                                  strat_1_rate_c_1 = log(2) / 6,
                                                  strat_1_rate_c_2 = log(2) / 6,
                                                  strat_1_rate_e_1 = log(2) / 6 * (2/3),
                                                  strat_1_rate_e_2 = log(2) / 6 * (2/3),
                                                  strat_0_rate_c_1 = log(2) / 10,
                                                  strat_0_rate_c_2 = log(2) / 10,
                                                  strat_0_rate_e_1 = log(2) / 10 * (2/3),
                                                  strat_0_rate_e_2 = log(2) / 10 * (2/3),
                                                  prop_strat_1 = 0.5,
                                                  label = "Moderate prog.",
                                                  effect = "PH",
                                                  scenario = "Scenario 1")


###### C. strong prognostic ##############

plot_homogeneous_strong_prognostic = plot_strat(t_end = 24, 
                                                t_star_1 = 6,
                                                t_star_0 = 6,
                                                strat_1_rate_c_1 = log(2) / 3,
                                                strat_1_rate_c_2 = log(2) / 3,
                                                strat_1_rate_e_1 = log(2) / 3 * (2/3),
                                                strat_1_rate_e_2 = log(2) / 3 * (2/3),
                                                strat_0_rate_c_1 = log(2) / 15,
                                                strat_0_rate_c_2 = log(2) / 15,
                                                strat_0_rate_e_1 = log(2) / 15 * (2/3),
                                                strat_0_rate_e_2 = log(2) / 15 * (2/3),
                                                prop_strat_1 = 0.5,
                                                label = "Strong prog.",
                                                effect = "PH",
                                                scenario = "Scenario 1")


##########################################
## 2. Quantitative difference (poor better)
##########################################

###### A. Non-prognostic #################

plot_quantitative_poor_better_non_prognostic = plot_strat(t_end = 24, 
                                                         t_star_1 = 6,
                                                         t_star_0 = 6,
                                                         strat_1_rate_c_1 = log(2) / 8,
                                                         strat_1_rate_c_2 = log(2) / 8,
                                                         strat_1_rate_e_1 = log(2) / 8 * 0.6,
                                                         strat_1_rate_e_2 = log(2) / 8 * 0.6,
                                                         strat_0_rate_c_1 = log(2) / 8,
                                                         strat_0_rate_c_2 = log(2) / 8,
                                                         strat_0_rate_e_1 = log(2) / 8 * (3/4),
                                                         strat_0_rate_e_2 = log(2) / 8 * (3/4),
                                                         prop_strat_1 = 0.5,
                                                         label = "Non-prog.",
                                                         effect = "PH",
                                                         scenario = "Scenario 2")


###### B. moderate prognostic ############

plot_quantitative_poor_better_moderate_prognostic = plot_strat(t_end = 24, 
                                                              t_star_1 = 6,
                                                              t_star_0 = 6,
                                                              strat_1_rate_c_1 = log(2) / 6,
                                                              strat_1_rate_c_2 = log(2) / 6,
                                                              strat_1_rate_e_1 = log(2) / 6 * 0.6,
                                                              strat_1_rate_e_2 = log(2) / 6 * 0.6,
                                                              strat_0_rate_c_1 = log(2) / 10,
                                                              strat_0_rate_c_2 = log(2) / 10,
                                                              strat_0_rate_e_1 = log(2) / 10 * (3/4),
                                                              strat_0_rate_e_2 = log(2) / 10 * (3/4),
                                                              prop_strat_1 = 0.5,
                                                              label = "Moderate prog.",
                                                              effect = "PH",
                                                              scenario = "Scenario 2")

###### C. strong prognostic ##############

plot_quantitative_poor_better_strong_prognostic = plot_strat(t_end = 24, 
                                                            t_star_1 = 6,
                                                            t_star_0 = 6,
                                                            strat_1_rate_c_1 = log(2) / 3,
                                                            strat_1_rate_c_2 = log(2) / 3,
                                                            strat_1_rate_e_1 = log(2) / 3 * 0.6,
                                                            strat_1_rate_e_2 = log(2) / 3 * 0.6,
                                                            strat_0_rate_c_1 = log(2) / 15,
                                                            strat_0_rate_c_2 = log(2) / 15,
                                                            strat_0_rate_e_1 = log(2) / 15 * (3/4),
                                                            strat_0_rate_e_2 = log(2) / 15 * (3/4),
                                                            prop_strat_1 = 0.5,
                                                            label = "Strong prog.",
                                                            effect = "PH",
                                                            scenario = "Scenario 2")


##########################################
## 3. Quantitative difference (poor worse)
##########################################

###### A. Non-prognostic #################

plot_quantitative_poor_worse_non_prognostic = plot_strat(t_end = 24, 
                                                        t_star_1 = 6,
                                                        t_star_0 = 6,
                                                        strat_1_rate_c_1 = log(2) / 8,
                                                        strat_1_rate_c_2 = log(2) / 8,
                                                        strat_1_rate_e_1 = log(2) / 8 * (3/4),
                                                        strat_1_rate_e_2 = log(2) / 8 * (3/4),
                                                        strat_0_rate_c_1 = log(2) / 8,
                                                        strat_0_rate_c_2 = log(2) / 8,
                                                        strat_0_rate_e_1 = log(2) / 8 * 0.6,
                                                        strat_0_rate_e_2 = log(2) / 8 * 0.6,
                                                        prop_strat_1 = 0.5,
                                                        label = "Non-prog.",
                                                        effect = "PH",
                                                        scenario = "Scenario 3")

###### B. moderate prognostic ############

plot_quantitative_poor_worse_moderate_prognostic = plot_strat(t_end = 24, 
                                                             t_star_1 = 6,
                                                             t_star_0 = 6,
                                                             strat_1_rate_c_1 = log(2) / 6,
                                                             strat_1_rate_c_2 = log(2) / 6,
                                                             strat_1_rate_e_1 = log(2) / 6 * (3/4),
                                                             strat_1_rate_e_2 = log(2) / 6 * (3/4),
                                                             strat_0_rate_c_1 = log(2) / 10,
                                                             strat_0_rate_c_2 = log(2) / 10,
                                                             strat_0_rate_e_1 = log(2) / 10 * 0.6,
                                                             strat_0_rate_e_2 = log(2) / 10 * 0.6,
                                                             prop_strat_1 = 0.5,
                                                             label = "Moderate prog.",
                                                             effect = "PH",
                                                             scenario = "Scenario 3")


###### C. strong prognostic ##############

plot_quantitative_poor_worse_strong_prognostic = plot_strat(t_end = 24, 
                                                           t_star_1 = 6,
                                                           t_star_0 = 6,
                                                           strat_1_rate_c_1 = log(2) / 3,
                                                           strat_1_rate_c_2 = log(2) / 3,
                                                           strat_1_rate_e_1 = log(2) / 3 * (3/4),
                                                           strat_1_rate_e_2 = log(2) / 3 * (3/4),
                                                           strat_0_rate_c_1 = log(2) / 15,
                                                           strat_0_rate_c_2 = log(2) / 15,
                                                           strat_0_rate_e_1 = log(2) / 15 * 0.6,
                                                           strat_0_rate_e_2 = log(2) / 15 * 0.6,
                                                           prop_strat_1 = 0.5,
                                                           label = "Strong prog.",
                                                           effect = "PH",
                                                           scenario = "Scenario 3")


##########################################
## 6. "Homogeneous (NPH)"
##########################################

###### A. Non-prognostic #################

plot_homogeneous_nph_non_prognostic = plot_strat(t_end = 24, 
                                                 t_star_1 = 6,
                                                 t_star_0 = 6,
                                                 strat_1_rate_c_1 = log(2) / 8,
                                                 strat_1_rate_c_2 = log(2) / 8,
                                                 strat_1_rate_e_1 = log(2) / 8 * (1),
                                                 strat_1_rate_e_2 = log(2) / 8 * (1/2),
                                                 strat_0_rate_c_1 = log(2) / 8,
                                                 strat_0_rate_c_2 = log(2) / 8,
                                                 strat_0_rate_e_1 = log(2) / 8 * 1,
                                                 strat_0_rate_e_2 = log(2) / 8 * (1/2),
                                                 prop_strat_1 = 0.5,
                                                 label = "Non-prog.",
                                                 effect = "NPH",
                                                 scenario = "Scenario 4")


###### B. moderate prognostic ############

plot_homogeneous_nph_moderate_prognostic = plot_strat(t_end = 24, 
                                                      t_star_1 = 6,
                                                      t_star_0 = 6,
                                                      strat_1_rate_c_1 = log(2) / 6,
                                                      strat_1_rate_c_2 = log(2) / 6,
                                                      strat_1_rate_e_1 = log(2) / 6 * (1),
                                                      strat_1_rate_e_2 = log(2) / 6 * (0.4),
                                                      strat_0_rate_c_1 = log(2) / 10,
                                                      strat_0_rate_c_2 = log(2) / 10,
                                                      strat_0_rate_e_1 = log(2) / 10 * (1),
                                                      strat_0_rate_e_2 = log(2) / 10 * (0.5),
                                                      prop_strat_1 = 0.5,
                                                      label = "Moderate prog.",
                                                      effect = "NPH",
                                                      scenario = "Scenario 4")


###### C. strong prognostic ##############

plot_homogeneous_nph_strong_prognostic = plot_strat(t_end = 24, 
                                                    t_star_1 = 6,
                                                    t_star_0 = 6,
                                                    strat_1_rate_c_1 = log(2) / 3,
                                                    strat_1_rate_c_2 = log(2) / 3,
                                                    strat_1_rate_e_1 = log(2) / 3 * (1),
                                                    strat_1_rate_e_2 = log(2) / 3 * (0.2),
                                                    strat_0_rate_c_1 = log(2) / 15,
                                                    strat_0_rate_c_2 = log(2) / 15,
                                                    strat_0_rate_e_1 = log(2) / 15 * (1),
                                                    strat_0_rate_e_2 = log(2) / 15 * (0.5),
                                                    prop_strat_1 = 0.5,
                                                    label = "Strong prog.",
                                                    effect = "NPH",
                                                    scenario = "Scenario 4")



###################################################
## 7. Quantitative difference (poor better) (NPH) 
###################################################

###### A. Non-prognostic #################

plot_quantitative_poor_better_nph_non_prognostic = plot_strat(t_end = 24, 
                                                              t_star_1 = 6,
                                                              t_star_0 = 6,
                                                              strat_1_rate_c_1 = log(2) / 8,
                                                              strat_1_rate_c_2 = log(2) / 8,
                                                              strat_1_rate_e_1 = log(2) / 8 * (1),
                                                              strat_1_rate_e_2 = log(2) / 8 * (0.3),
                                                              strat_0_rate_c_1 = log(2) / 8,
                                                              strat_0_rate_c_2 = log(2) / 8,
                                                              strat_0_rate_e_1 = log(2) / 8 * 1,
                                                              strat_0_rate_e_2 = log(2) / 8 * (0.7),
                                                              prop_strat_1 = 0.5,
                                                              label = "Non-prog.",
                                                              effect = "NPH",
                                                              scenario = "Scenario 5")


###### B. moderate prognostic ############

plot_quantitative_poor_better_nph_moderate_prognostic = plot_strat(t_end = 24, 
                                                                   t_star_1 = 6,
                                                                   t_star_0 = 6,
                                                                   strat_1_rate_c_1 = log(2) / 6,
                                                                   strat_1_rate_c_2 = log(2) / 6,
                                                                   strat_1_rate_e_1 = log(2) / 6 * (1),
                                                                   strat_1_rate_e_2 = log(2) / 6 * (0.3),
                                                                   strat_0_rate_c_1 = log(2) / 10,
                                                                   strat_0_rate_c_2 = log(2) / 10,
                                                                   strat_0_rate_e_1 = log(2) / 10 * (1),
                                                                   strat_0_rate_e_2 = log(2) / 10 * (0.7),
                                                                   prop_strat_1 = 0.5,
                                                                   label = "Moderate prog.",
                                                                   effect = "NPH",
                                                                   scenario = "Scenario 5")


###### C. strong prognostic ##############

plot_quantitative_poor_better_nph_strong_prognostic = plot_strat(t_end = 24, 
                                                                 t_star_1 = 6,
                                                                 t_star_0 = 6,
                                                                 strat_1_rate_c_1 = log(2) / 3,
                                                                 strat_1_rate_c_2 = log(2) / 3,
                                                                 strat_1_rate_e_1 = log(2) / 3 * (1),
                                                                 strat_1_rate_e_2 = log(2) / 3 * (0.1),
                                                                 strat_0_rate_c_1 = log(2) / 15,
                                                                 strat_0_rate_c_2 = log(2) / 15,
                                                                 strat_0_rate_e_1 = log(2) / 15 * (1),
                                                                 strat_0_rate_e_2 = log(2) / 15 * (0.8),
                                                                 prop_strat_1 = 0.5,
                                                                 label = "Strong prog.",
                                                                 effect = "NPH",
                                                                 scenario = "Scenario 5")


###################################################
## 8. Quantitative difference (poor worse)  (NPH)
###################################################

###### A. Non-prognostic #################

plot_quantitative_poor_worse_nph_non_prognostic = plot_strat(t_end = 24, 
                                                             t_star_1 = 6,
                                                             t_star_0 = 6,
                                                             strat_1_rate_c_1 = log(2) / 8,
                                                             strat_1_rate_c_2 = log(2) / 8,
                                                             strat_1_rate_e_1 = log(2) / 8 * (1),
                                                             strat_1_rate_e_2 = log(2) / 8 * (0.7),
                                                             strat_0_rate_c_1 = log(2) / 8,
                                                             strat_0_rate_c_2 = log(2) / 8,
                                                             strat_0_rate_e_1 = log(2) / 8 * 1,
                                                             strat_0_rate_e_2 = log(2) / 8 * (0.3),
                                                             prop_strat_1 = 0.5,
                                                             label = "Non-prog.",
                                                             effect = "NPH",
                                                             scenario = "Scenario 6")


###### B. moderate prognostic ############

plot_quantitative_poor_worse_nph_moderate_prognostic = plot_strat(t_end = 24, 
                                                                  t_star_1 = 6,
                                                                  t_star_0 = 6,
                                                                  strat_1_rate_c_1 = log(2) / 6,
                                                                  strat_1_rate_c_2 = log(2) / 6,
                                                                  strat_1_rate_e_1 = log(2) / 6 * (1),
                                                                  strat_1_rate_e_2 = log(2) / 6 * (0.7),
                                                                  strat_0_rate_c_1 = log(2) / 10,
                                                                  strat_0_rate_c_2 = log(2) / 10,
                                                                  strat_0_rate_e_1 = log(2) / 10 * (1),
                                                                  strat_0_rate_e_2 = log(2) / 10 * (0.3),
                                                                  prop_strat_1 = 0.5,
                                                                  label = "Moderate prog.",
                                                                  effect = "NPH",
                                                                  scenario = "Scenario 6")


###### C. strong prognostic ##############

plot_quantitative_poor_worse_nph_strong_prognostic = plot_strat(t_end = 24, 
                                                                t_star_1 = 6,
                                                                t_star_0 = 6,
                                                                strat_1_rate_c_1 = log(2) / 3,
                                                                strat_1_rate_c_2 = log(2) / 3,
                                                                strat_1_rate_e_1 = log(2) / 3 * (1),
                                                                strat_1_rate_e_2 = log(2) / 3 * (0.7),
                                                                strat_0_rate_c_1 = log(2) / 15,
                                                                strat_0_rate_c_2 = log(2) / 15,
                                                                strat_0_rate_e_1 = log(2) / 15 * (1),
                                                                strat_0_rate_e_2 = log(2) / 15 * (0.3),
                                                                prop_strat_1 = 0.5,
                                                                label = "Strong prog.",
                                                                effect = "NPH",
                                                                scenario = "Scenario 6")




##########################################
## 0. Null effect
##########################################

###### A. Non-prognostic #################
plot_null_effect_non_prognostic = plot_strat(t_end = 24, 
                                             t_star_1 = 6,
                                             t_star_0 = 6,
                                             strat_1_rate_c_1 = log(2) / 8,
                                             strat_1_rate_c_2 = log(2) / 8,
                                             strat_1_rate_e_1 = log(2) / 8 * (1),
                                             strat_1_rate_e_2 = log(2) / 8 * (1),
                                             strat_0_rate_c_1 = log(2) / 8,
                                             strat_0_rate_c_2 = log(2) / 8,
                                             strat_0_rate_e_1 = log(2) / 8 * (1),
                                             strat_0_rate_e_2 = log(2) / 8 * (1),
                                             prop_strat_1 = 0.5,
                                             label = "Non-prog.",
                                             effect = "Null effect",
                                             scenario = "Scenario 7")

###### B. moderate prognostic ############
plot_null_effect_moderate_prognostic = plot_strat(t_end = 24, 
                                                  t_star_1 = 6,
                                                  t_star_0 = 6,
                                                  strat_1_rate_c_1 = log(2) / 6,
                                                  strat_1_rate_c_2 = log(2) / 6,
                                                  strat_1_rate_e_1 = log(2) / 6 * (1),
                                                  strat_1_rate_e_2 = log(2) / 6 * (1),
                                                  strat_0_rate_c_1 = log(2) / 10,
                                                  strat_0_rate_c_2 = log(2) / 10,
                                                  strat_0_rate_e_1 = log(2) / 10 * (1),
                                                  strat_0_rate_e_2 = log(2) / 10 * (1),
                                                  prop_strat_1 = 0.5,
                                                  label = "Moderate prog.",
                                                  effect = "Null effect",
                                                  scenario = "Scenario 7")

###### C. strong prognostic ##############
plot_null_effect_strong_prognostic = plot_strat(t_end = 24, 
                                                t_star_1 = 6,
                                                t_star_0 = 6,
                                                strat_1_rate_c_1 = log(2) / 3,
                                                strat_1_rate_c_2 = log(2) / 3,
                                                strat_1_rate_e_1 = log(2) / 3 * (1),
                                                strat_1_rate_e_2 = log(2) / 3 * (1),
                                                strat_0_rate_c_1 = log(2) / 15,
                                                strat_0_rate_c_2 = log(2) / 15,
                                                strat_0_rate_e_1 = log(2) / 15 * (1),
                                                strat_0_rate_e_2 = log(2) / 15 * (1),
                                                prop_strat_1 = 0.5,
                                                label = "Strong prog.",
                                                effect = "Null effect",
                                                scenario = "Scenario 7")



##########################################
## 4. Qualitative difference (poor better)
##########################################

###### A. Non-prognostic #################

plot_qualitative_poor_better_non_prognostic = plot_strat(t_end = 24, 
                                                         t_star_1 = 11,
                                                         t_star_0 = 11,
                                                         strat_1_rate_c_1 = log(2) / 8,
                                                         strat_1_rate_c_2 = log(2) / 8,
                                                         strat_1_rate_e_1 = log(2) / 8 * (2/3),
                                                         strat_1_rate_e_2 = log(2) / 8 * 0.72,
                                                         strat_0_rate_c_1 = log(2) / 8,
                                                         strat_0_rate_c_2 = log(2) / 8,
                                                         strat_0_rate_e_1 = log(2) / 8 * 1.5,
                                                         strat_0_rate_e_2 = log(2) / 8 * 2.2,
                                                         prop_strat_1 = 0.5,
                                                         label = "Non-prog.",
                                                         effect = "Null effect",
                                                         scenario = "Scenario 8")


###### B. moderate prognostic ############

plot_qualitative_poor_better_moderate_prognostic = plot_strat(t_end = 24, 
                                                              t_star_1 = 11,
                                                              t_star_0 = 11,
                                                              strat_1_rate_c_1 = log(2) / 6,
                                                              strat_1_rate_c_2 = log(2) / 6,
                                                              strat_1_rate_e_1 = log(2) / 6 * (2/3),
                                                              strat_1_rate_e_2 = log(2) / 6 * (2/3),
                                                              strat_0_rate_c_1 = log(2) / 10,
                                                              strat_0_rate_c_2 = log(2) / 10,
                                                              strat_0_rate_e_1 = log(2) / 10 * 1.5,
                                                              strat_0_rate_e_2 = log(2) / 10 * 1.5,
                                                              prop_strat_1 = 0.5,
                                                              label = "Moderate prog.",
                                                              effect = "Null effect",
                                                              scenario = "Scenario 8")


###### C. strong prognostic ##############

plot_qualitative_poor_better_strong_prognostic = plot_strat(t_end = 24, 
                                                            t_star_1 = 6,
                                                            t_star_0 = 8,
                                                            strat_1_rate_c_1 = log(2) / 3,
                                                            strat_1_rate_c_2 = log(2) / 3,
                                                            strat_1_rate_e_1 = log(2) / 3 * (2/3),
                                                            strat_1_rate_e_2 = log(2) / 3 * 0.3,
                                                            strat_0_rate_c_1 = log(2) / 15,
                                                            strat_0_rate_c_2 = log(2) / 15,
                                                            strat_0_rate_e_1 = log(2) / 15 * 2.1,
                                                            strat_0_rate_e_2 = log(2) / 15 * 1.2,
                                                            prop_strat_1 = 0.5,
                                                            label = "Strong prog.",
                                                            effect = "Null effect",
                                                            scenario = "Scenario 8")




##########################################
## 5. Qualitative difference (poor worse)
##########################################

###### A. Non-prognostic #################

plot_qualitative_poor_worse_non_prognostic = plot_strat(t_end = 24, 
                                                        t_star_1 = 9,
                                                        t_star_0 = 9,
                                                        strat_1_rate_c_1 = log(2) / 8,
                                                        strat_1_rate_c_2 = log(2) / 8,
                                                        strat_1_rate_e_1 = log(2) / 8 * 1.8,
                                                        strat_1_rate_e_2 = log(2) / 8 * 2,
                                                        strat_0_rate_c_1 = log(2) / 8,
                                                        strat_0_rate_c_2 = log(2) / 8,
                                                        strat_0_rate_e_1 = log(2) / 8 * 0.5,
                                                        strat_0_rate_e_2 = log(2) / 8 * 0.8,
                                                        prop_strat_1 = 0.5,
                                                        label = "Non-prog.",
                                                        effect = "Null effect",
                                                        scenario = "Scenario 9")


###### B. moderate prognostic ############

plot_qualitative_poor_worse_moderate_prognostic = plot_strat(t_end = 24, 
                                                             t_star_1 = 11,
                                                             t_star_0 = 11,
                                                             strat_1_rate_c_1 = log(2) / 6,
                                                             strat_1_rate_c_2 = log(2) / 6,
                                                             strat_1_rate_e_1 = log(2) / 6 * 1.5,
                                                             strat_1_rate_e_2 = log(2) / 6 * 2.4,
                                                             strat_0_rate_c_1 = log(2) / 10,
                                                             strat_0_rate_c_2 = log(2) / 10,
                                                             strat_0_rate_e_1 = log(2) / 10 * 0.65,
                                                             strat_0_rate_e_2 = log(2) / 10 * 0.95,
                                                             prop_strat_1 = 0.5,
                                                             label = "Moderate prog.",
                                                             effect = "Null effect",
                                                             scenario = "Scenario 9")


###### C. strong prognostic ##############

plot_qualitative_poor_worse_strong_prognostic = plot_strat(t_end = 24, 
                                                           t_star_1 = 6,
                                                           t_star_0 = 6,
                                                           strat_1_rate_c_1 = log(2) / 3,
                                                           strat_1_rate_c_2 = log(2) / 3,
                                                           strat_1_rate_e_1 = log(2) / 3 * 1.4,
                                                           strat_1_rate_e_2 = log(2) / 3 * (2),
                                                           strat_0_rate_c_1 = log(2) / 15,
                                                           strat_0_rate_c_2 = log(2) / 15,
                                                           strat_0_rate_e_1 = log(2) / 15 * 0.675,
                                                           strat_0_rate_e_2 = log(2) / 15 * 1.075,
                                                           prop_strat_1 = 0.5,
                                                           label = "Strong prog.",
                                                           effect = "Null effect",
                                                           scenario = "Scenario 9")



# --- PLOTS --- #

#PH

data_scenario_PH = rbind(plot_homogeneous_non_prognostic,
                         plot_homogeneous_moderate_prognostic,
                         plot_homogeneous_strong_prognostic,
                         plot_quantitative_poor_better_non_prognostic,
                         plot_quantitative_poor_better_moderate_prognostic,
                         plot_quantitative_poor_better_strong_prognostic,
                         plot_quantitative_poor_worse_non_prognostic,
                         plot_quantitative_poor_worse_moderate_prognostic,
                         plot_quantitative_poor_worse_strong_prognostic)

data_scenario_PH$prognosis = factor(data_scenario_PH$prognosis, levels = c("Non-prog.","Moderate prog.","Strong prog."))


plot_scenarios_PH = 
  ggplot(data_scenario_PH, aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank())

ggsave("plot_scenarios_PH.pdf", plot_scenarios_PH, width = 8, height = 12, dpi = 300)


#NPH


data_scenario_NPH = rbind(plot_homogeneous_nph_non_prognostic,
                          plot_homogeneous_nph_moderate_prognostic,
                          plot_homogeneous_nph_strong_prognostic,
                          plot_quantitative_poor_better_nph_non_prognostic,
                          plot_quantitative_poor_better_nph_moderate_prognostic,
                          plot_quantitative_poor_better_nph_strong_prognostic,
                          plot_quantitative_poor_worse_nph_non_prognostic,
                          plot_quantitative_poor_worse_nph_moderate_prognostic,
                          plot_quantitative_poor_worse_nph_strong_prognostic)

data_scenario_NPH$prognosis = factor(data_scenario_NPH$prognosis, levels = c("Non-prog.","Moderate prog.","Strong prog."))


plot_scenarios_NPH = 
  ggplot(data_scenario_NPH, aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank())

ggsave("plot_scenarios_NPH.pdf", plot_scenarios_NPH, width = 8, height = 12, dpi = 300)


#NULL EFFECT

data_scenario_null = rbind(plot_null_effect_non_prognostic,
                           plot_null_effect_moderate_prognostic,
                           plot_null_effect_strong_prognostic,
                           plot_qualitative_poor_better_non_prognostic,
                           plot_qualitative_poor_better_moderate_prognostic,
                           plot_qualitative_poor_better_strong_prognostic,
                           plot_qualitative_poor_worse_non_prognostic,
                           plot_qualitative_poor_worse_moderate_prognostic,
                           plot_qualitative_poor_worse_strong_prognostic)

data_scenario_null$prognosis = factor(data_scenario_null$prognosis, levels = c("Non-prog.","Moderate prog.","Strong prog."))


plot_scenarios_null = 
  ggplot(data_scenario_null, aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank())

ggsave("plot_scenarios_null.pdf", plot_scenarios_null, width = 8, height = 12, dpi = 300)
