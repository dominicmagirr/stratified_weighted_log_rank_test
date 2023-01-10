library(dplyr)
library(ggplot2)
library(survival)

#This file produces Figures 6, 7 and 8.

#############################################
#Set the output, data and code folders paths#
#############################################

data_path = "./Data/"
output_path = "./Outputs/"
code_path = "./Code/"

source(paste0(code_path,"scenarios_list.R"))



#Function needed to plot different stratified scenarios
plot_strat <- function(scenario, t_end = 24, prop_strat_1 = 0.5){
  
  
  #Function needed to compute RMST
  surv_2 <- function(t, lambdas, changepoint){
    
    p_1 <- exp(-lambdas[1] * t)
    p_2 <- exp(-lambdas[1] * changepoint) * exp(-lambdas[2] * pmax(0, (t - changepoint)))
    
    (t < changepoint) * p_1 + (t >= changepoint) * p_2
    
  }
  
  #Function needed to compute RMST
  get_RMST = function(lambdas, changepoint, tau){
    
    
    rmst = integrate(surv_2, 
                     lower = 0, 
                     upper = tau,
                     lambdas = lambdas,
                     changepoint = changepoint)$value
    
    return(rmst)
    
  }
  
  t_seq <- seq(0, t_end, length.out = 100)
  t_length <- length(t_seq)
  
  s_c_1 <- surv_2(t_seq, lambdas = scenario$strata_1_c, changepoint = scenario$changepoint_1)
  s_e_1 <- surv_2(t_seq, lambdas = scenario$strata_1_e, changepoint = scenario$changepoint_1)
  s_c_2 <- surv_2(t_seq, lambdas = scenario$strata_2_c, changepoint = scenario$changepoint_2)
  s_e_2 <- surv_2(t_seq, lambdas = scenario$strata_2_e, changepoint = scenario$changepoint_2)
  
  s_c <- s_c_1 * prop_strat_1 + s_c_2 * (1 - prop_strat_1)
  s_e <- s_e_1 * prop_strat_1 + s_e_2 * (1 - prop_strat_1)
  
  rmst_c_1 <- get_RMST(scenario$strata_1_c, scenario$changepoint_1, t_end)
  rmst_e_1 <- get_RMST(scenario$strata_1_e, scenario$changepoint_1, t_end)
  rmst_c_2 <- get_RMST(scenario$strata_2_c, scenario$changepoint_2, t_end)
  rmst_e_2 <- get_RMST(scenario$strata_2_e, scenario$changepoint_2, t_end)
  
  rmst_c <- rmst_c_1 * prop_strat_1 + rmst_c_2 * (1 - prop_strat_1) 
  rmst_e <- rmst_e_1 * prop_strat_1 + rmst_e_2 * (1 - prop_strat_1) 
  
  df_1 <- data.frame(t_seq = rep(t_seq, 2), 
                     surv = c(s_c_1, s_e_1),
                     arm = rep(c("Control", "Experimental"), each = t_length),
                     stratum = "First stratum",
                     rmst = rep(c(rmst_c_1, rmst_e_1) %>% round(1), each = t_length))
  
  df_2 <- data.frame(t_seq = rep(t_seq, 2), 
                     surv = c(s_c_2, s_e_2),
                     arm = rep(c("Control", "Experimental"), each = t_length),
                     stratum = "Second stratum",
                     rmst = rep(c(rmst_c_2, rmst_e_2) %>% round(1), each = t_length))
   
  df_m <- data.frame(t_seq = rep(t_seq, 2), 
                     surv = c(s_c, s_e),
                     arm = rep(c("Control", "Experimental"), each = t_length),
                     stratum = "Overall",
                     rmst = rep(c(rmst_c, rmst_e) %>% round(1), each = t_length))
  
  
  plots_df = rbind(df_1, df_2, df_m) %>%
    mutate(prognosis = scenario$label,
           effect = scenario$effect,
           scenario = scenario$scenario_label)
  
  return(plots_df)
  

}


################################
## Proportional Hazards Scenarios
################################
plot_1_A <- plot_strat(scenario_1_A)
plot_1_B <- plot_strat(scenario_1_B)
plot_1_C <- plot_strat(scenario_1_C)
##
plot_2_A <- plot_strat(scenario_2_A)
plot_2_B <- plot_strat(scenario_2_B)
plot_2_C <- plot_strat(scenario_2_C)
##
plot_3_A <- plot_strat(scenario_3_A)
plot_3_B <- plot_strat(scenario_3_B)
plot_3_C <- plot_strat(scenario_3_C)
################################
scenario_PH = rbind(plot_1_A,
                    plot_1_B,
                    plot_1_C,
                    ##
                    plot_2_A,
                    plot_2_B,
                    plot_2_C,
                    ##
                    plot_3_A,
                    plot_3_B,
                    plot_3_C)


scenario_PH$prognosis = factor(scenario_PH$prognosis,
                               levels = c("Non-prog.","Moderate prog.","Strong prog."))

scenario_PH$stratum = factor(scenario_PH$stratum,
                               levels = c("First stratum", "Second stratum", "Overall"))

plot_scenarios_PH <- ggplot(scenario_PH,
       aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank()) +
  geom_text(scenario_PH %>% filter(arm == "Control") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.82, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE) +
  geom_text(scenario_PH %>% filter(arm == "Experimental") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.95, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE)


ggsave(paste0(output_path,"Figure_6.pdf"), plot_scenarios_PH, width = 8, height = 12, dpi = 300)


#####################################
## Non-Proportional Hazards Scenarios
#####################################
plot_4_A <- plot_strat(scenario_4_A)
plot_4_B <- plot_strat(scenario_4_B)
plot_4_C <- plot_strat(scenario_4_C)
##
plot_5_A <- plot_strat(scenario_5_A)
plot_5_B <- plot_strat(scenario_5_B)
plot_5_C <- plot_strat(scenario_5_C)
##
plot_6_A <- plot_strat(scenario_6_A)
plot_6_B <- plot_strat(scenario_6_B)
plot_6_C <- plot_strat(scenario_6_C)
################################
scenario_NPH = rbind(plot_4_A,
                     plot_4_B,
                     plot_4_C,
                     ##
                     plot_5_A,
                     plot_5_B,
                     plot_5_C,
                     ##
                     plot_6_A,
                     plot_6_B,
                     plot_6_C)


scenario_NPH$prognosis = factor(scenario_NPH$prognosis,
                                levels = c("Non-prog.","Moderate prog.","Strong prog."))

scenario_NPH$stratum = factor(scenario_NPH$stratum,
                              levels = c("First stratum", "Second stratum", "Overall"))

plot_scenarios_NPH <- ggplot(scenario_NPH,
                             aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank()) +
  geom_text(scenario_NPH %>% filter(arm == "Control") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.82, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE) +
  geom_text(scenario_NPH %>% filter(arm == "Experimental") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.95, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE)


ggsave(paste0(output_path,"Figure_7.pdf"), plot_scenarios_NPH, width = 8, height = 12, dpi = 300)


#####################################
## Scenarios with null overall effect
#####################################
plot_7_A <- plot_strat(scenario_7_A)
plot_7_B <- plot_strat(scenario_7_B)
plot_7_C <- plot_strat(scenario_7_C)
##
plot_8_A <- plot_strat(scenario_8_A)
plot_8_B <- plot_strat(scenario_8_B)
plot_8_C <- plot_strat(scenario_8_C)
##
plot_9_A <- plot_strat(scenario_9_A)
plot_9_B <- plot_strat(scenario_9_B)
plot_9_C <- plot_strat(scenario_9_C)
################################
scenario_null = rbind(plot_7_A,
                      plot_7_B,
                      plot_7_C,
                      ##
                      plot_8_A,
                      plot_8_B,
                      plot_8_C,
                      ##
                      plot_9_A,
                      plot_9_B,
                      plot_9_C)


scenario_null$prognosis = factor(scenario_null$prognosis,
                                 levels = c("Non-prog.","Moderate prog.","Strong prog."))

scenario_null$stratum = factor(scenario_null$stratum,
                               levels = c("First stratum", "Second stratum", "Overall"))

plot_scenarios_null <- ggplot(scenario_null,
                              aes(x = t_seq, y = surv, color = arm)) +
  geom_line(size = 1) +
  scale_y_continuous("Survival", limits = c(0,1)) +
  scale_x_continuous("Time", limits = c(0,25)) +
  theme_bw() +
  facet_grid(scenario + prognosis ~ stratum ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank()) +
  geom_text(scenario_null %>% filter(arm == "Control") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.82, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE) +
  geom_text(scenario_null %>% filter(arm == "Experimental") %>% group_by(scenario, prognosis, stratum) %>% slice(1), 
            mapping = aes(x=20, y=0.95, label=paste("RMST = ",rmst)), size=3.25, show.legend = FALSE)


ggsave(paste0(output_path,"Figure_8.pdf"), plot_scenarios_null, width = 8, height = 12, dpi = 300)








