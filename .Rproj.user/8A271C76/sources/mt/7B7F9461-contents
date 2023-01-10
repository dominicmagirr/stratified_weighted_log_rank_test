library(tidyr)
##############################################
# load nphRCT (https://github.com/cran/nphRCT)
devtools::load_all("~/nphRCT/")

source("scenarios_list.R")
recruitment_model <- list(rec_model="power", rec_period = 9, rec_power = 1)

sim_one_trial <- function(t_end = 24, scenario, t_star = 12, s_star = NULL){
  
  ################################################
  ## strata 1 data
  sim_data_1 <- nphRCT::sim_events_delay(
    event_model=list(
      duration_c = c(scenario$changepoint_1, t_end),
      duration_e = c(scenario$changepoint_1, t_end),
      lambda_c = scenario$strata_1_c,
      lambda_e = scenario$strata_1_e
    ),
    recruitment_model=recruitment_model,
    n_c=86,
    n_e=86,
    max_cal_t = t_end
  ) %>% 
    mutate(strata = "first")
  ################################################
  ## strata 2 data
  sim_data_2 <- nphRCT::sim_events_delay(
    event_model=list(
      duration_c = c(scenario$changepoint_2, t_end),
      duration_e = c(scenario$changepoint_2, t_end),
      lambda_c = scenario$strata_2_c,
      lambda_e = scenario$strata_2_e
    ),
    recruitment_model=recruitment_model,
    n_c=86,
    n_e=86,
    max_cal_t = t_end
  ) %>% 
    mutate(strata = "second")
  ################################################
  ## full data
  sim_data <- rbind(sim_data_1, sim_data_2)
  
  ###############################################
  ## wlrt
  res_strata_0 <- nphRCT::wlrt(Surv(event_time, event_status) ~ group + strata(strata), data = sim_data, method = "mw", t_star = 0)
  res_strata <- nphRCT::wlrt(Surv(event_time, event_status) ~ group + strata(strata), data = sim_data, method = "mw", t_star = t_star, s_star = s_star)
  res_unstrata_0 <- nphRCT::wlrt(Surv(event_time, event_status) ~ group, data = sim_data, method = "mw", t_star = 0)
  res_unstrata <- nphRCT::wlrt(Surv(event_time, event_status) ~ group, data = sim_data, method = "mw", t_star = t_star, s_star = s_star)
  
  #### main Z tests
  z <- res_unstrata_0$z
  z_w <- res_unstrata$z
  tilde_z <- res_strata_0$combined$z
  tilde_z_w_z <- res_strata$combined$z
  
  #### create tilde_z_w_u
  tilde_z_w_u <- sum(res_strata$by_strata$u) / sqrt(sum(res_strata$by_strata$v_u))
  
  #### create tilde_z_w_n
  n <- c(172, 172)
  v_u_w <- res_strata$by_strata$v_u
  tilde_z_w_n <- sum(res_strata$by_strata$u / v_u_w * n) / sqrt(sum(n^2 / v_u_w))
  
  
  #### create tilde_z_n
  v_u <- res_strata_0$by_strata$v_u
  tilde_z_n <- sum(res_strata_0$by_strata$u / v_u * n) / sqrt(sum(n^2 / v_u))
  
  
  #### return z's
  data.frame(z = z,
             z_w = z_w,
             tilde_z = tilde_z,
             tilde_z_w_u = tilde_z_w_u,
             tilde_z_w_z = tilde_z_w_z,
             tilde_z_n = tilde_z_n,
             tilde_z_w_n = tilde_z_w_n)
  
}



#################################
#Get empirical power
get_power <- function(n_sim = 100, t_end = 24, scenario, t_star = 12, s_star = NULL){
  
  sim_res <- purrr::map_df(rep(t_end, n_sim), sim_one_trial, scenario = scenario, t_star = t_star, s_star = s_star)
  apply(sim_res, 2, function(x) mean(x < qnorm(0.025)))
  
  
}

###################
set.seed(124543)
###################
power_1_A <- get_power(n_sim = 10000, scenario = scenario_1_A)
power_1_B <- get_power(n_sim = 10000, scenario = scenario_1_B)
power_1_C <- get_power(n_sim = 10000, scenario = scenario_1_C)

power_2_A <- get_power(n_sim = 10000, scenario = scenario_2_A)
power_2_B <- get_power(n_sim = 10000, scenario = scenario_2_B)
power_2_C <- get_power(n_sim = 10000, scenario = scenario_2_C)

power_3_A <- get_power(n_sim = 10000, scenario = scenario_3_A)
power_3_B <- get_power(n_sim = 10000, scenario = scenario_3_B)
power_3_C <- get_power(n_sim = 10000, scenario = scenario_3_C)

power_4_A <- get_power(n_sim = 10000, scenario = scenario_4_A)
power_4_B <- get_power(n_sim = 10000, scenario = scenario_4_B)
power_4_C <- get_power(n_sim = 10000, scenario = scenario_4_C)

power_5_A <- get_power(n_sim = 10000, scenario = scenario_5_A)
power_5_B <- get_power(n_sim = 10000, scenario = scenario_5_B)
power_5_C <- get_power(n_sim = 10000, scenario = scenario_5_C)

power_6_A <- get_power(n_sim = 10000, scenario = scenario_6_A)
power_6_B <- get_power(n_sim = 10000, scenario = scenario_6_B)
power_6_C <- get_power(n_sim = 10000, scenario = scenario_6_C)

power_7_A <- get_power(n_sim = 10000, scenario = scenario_7_A)
power_7_B <- get_power(n_sim = 10000, scenario = scenario_7_B)
power_7_C <- get_power(n_sim = 10000, scenario = scenario_7_C)

power_8_A <- get_power(n_sim = 10000, scenario = scenario_8_A)
power_8_B <- get_power(n_sim = 10000, scenario = scenario_8_B)
power_8_C <- get_power(n_sim = 10000, scenario = scenario_8_C)

power_9_A <- get_power(n_sim = 10000, scenario = scenario_9_A)
power_9_B <- get_power(n_sim = 10000, scenario = scenario_9_B)
power_9_C <- get_power(n_sim = 10000, scenario = scenario_9_C)


################################################################

data_1 <- rbind(power_1_A, power_1_B, power_1_C) %>% 
  cbind(treatment_effect = rep("Scenario 1 (Homogeneous HRs)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_2 <- rbind(power_2_A, power_2_B, power_2_C) %>% 
  cbind(treatment_effect = rep("Scenario 2 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_3 <- rbind(power_3_A, power_3_B, power_3_C) %>% 
  cbind(treatment_effect = rep("Scenario 3 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")


################################################################

data_4 <- rbind(power_4_A, power_4_B, power_4_C) %>% 
  cbind(treatment_effect = rep("Scenario 4 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_5 <- rbind(power_5_A, power_5_B, power_5_C) %>% 
  cbind(treatment_effect = rep("Scenario 5 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")


data_6 <- rbind(power_6_A, power_6_B, power_6_C) %>% 
  cbind(treatment_effect = rep("Scenario 6 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

###################################################################

data_7 <- rbind(power_7_A, power_7_B, power_7_C) %>% 
  cbind(treatment_effect = rep("Scenario 7 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_8 <- rbind(power_8_A, power_8_B, power_8_C) %>% 
  cbind(treatment_effect = rep("Scenario 8 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_9 <- rbind(power_9_A, power_9_B, power_9_C) %>% 
  cbind(treatment_effect = rep("Scenario 9 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

###################################################################

null_ph_nph <- rbind(cbind(rbind(data_1, data_2, data_3),
                           pattern = "Proportional Hazards"),
                     cbind(rbind(data_4, data_5, data_6),
                           pattern = "Anticipated Delayed Effects"),
                     cbind(rbind(data_7, data_8, data_9),
                           pattern = "Null Overall Treatment Effect")) 

null_ph_nph$prognostic_effect <- factor(null_ph_nph$prognostic_effect,
                                        levels = c("None", "Moderate", "Strong"))

null_ph_nph$pattern <- factor(null_ph_nph$pattern,
                              levels = c("Proportional Hazards", "Anticipated Delayed Effects", "Null Overall Treatment Effect"))

save(null_ph_nph, file = "null_ph_nph.RData")

#Plos with all tests
plot_power = ggplot(data = null_ph_nph,
                    mapping = aes(x = prognostic_effect,
                                  y = as.numeric(as.character(value)),
                                  fill = test)) +
  facet_wrap(pattern ~ treatment_effect, scales = "free_y") +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_brewer(palette="Accent",
                    labels = c(expression("Strat. LR " (tilde(Z))), 
                               expression("Strat. LR from Mehrotra et al. " (Z^n)),
                               expression("Strat. WLR on sample size scale " (tilde(Z)^wn)),
                               expression("Strat. WLR on U-statistic scale " (tilde(Z)^wu)),
                               expression("Strat. WLR on Z-statistic scale " (tilde(Z)^wz)),
                               expression("Unstrat. LR " (Z)),
                               expression("Unstrat. WLR " (Z^w)))) +
  theme_bw() +
  #geom_hline(yintercept = 0.9, linetype = 2) + 
  geom_hline(yintercept = 0.025, linetype = 2) +
  scale_y_continuous("Power") +
  xlab("Prognostic effect") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text.align = 0)


plot_power
# Labels of the Z statistics with equations
# scale_fill_discrete(labels = c(expression(tilde(Z)^LR), 
#                                expression("Z"^n),
#                                expression(tilde(Z)^wn),
#                                expression(tilde(Z)^wu),
#                                expression(tilde(Z)^wz),
#                                expression(Z^LR),
#                                expression("Z"^w)))


ggsave("plot_tests_2.pdf", plot_power, width = 10, height = 7, dpi = 300)



#####################################################
## Repeat for s_star = 0.5
#####################################################


###################
set.seed(124543)
###################
power_1_A <- get_power(n_sim = 10000, scenario = scenario_1_A, t_star = NULL, s_star = 0.5)
power_1_B <- get_power(n_sim = 10000, scenario = scenario_1_B, t_star = NULL, s_star = 0.5)
power_1_C <- get_power(n_sim = 10000, scenario = scenario_1_C, t_star = NULL, s_star = 0.5)

power_2_A <- get_power(n_sim = 10000, scenario = scenario_2_A, t_star = NULL, s_star = 0.5)
power_2_B <- get_power(n_sim = 10000, scenario = scenario_2_B, t_star = NULL, s_star = 0.5)
power_2_C <- get_power(n_sim = 10000, scenario = scenario_2_C, t_star = NULL, s_star = 0.5)

power_3_A <- get_power(n_sim = 10000, scenario = scenario_3_A, t_star = NULL, s_star = 0.5)
power_3_B <- get_power(n_sim = 10000, scenario = scenario_3_B, t_star = NULL, s_star = 0.5)
power_3_C <- get_power(n_sim = 10000, scenario = scenario_3_C, t_star = NULL, s_star = 0.5)

power_4_A <- get_power(n_sim = 10000, scenario = scenario_4_A, t_star = NULL, s_star = 0.5)
power_4_B <- get_power(n_sim = 10000, scenario = scenario_4_B, t_star = NULL, s_star = 0.5)
power_4_C <- get_power(n_sim = 10000, scenario = scenario_4_C, t_star = NULL, s_star = 0.5)

power_5_A <- get_power(n_sim = 10000, scenario = scenario_5_A, t_star = NULL, s_star = 0.5)
power_5_B <- get_power(n_sim = 10000, scenario = scenario_5_B, t_star = NULL, s_star = 0.5)
power_5_C <- get_power(n_sim = 10000, scenario = scenario_5_C, t_star = NULL, s_star = 0.5)

power_6_A <- get_power(n_sim = 10000, scenario = scenario_6_A, t_star = NULL, s_star = 0.5)
power_6_B <- get_power(n_sim = 10000, scenario = scenario_6_B, t_star = NULL, s_star = 0.5)
power_6_C <- get_power(n_sim = 10000, scenario = scenario_6_C, t_star = NULL, s_star = 0.5)

power_7_A <- get_power(n_sim = 10000, scenario = scenario_7_A, t_star = NULL, s_star = 0.5)
power_7_B <- get_power(n_sim = 10000, scenario = scenario_7_B, t_star = NULL, s_star = 0.5)
power_7_C <- get_power(n_sim = 10000, scenario = scenario_7_C, t_star = NULL, s_star = 0.5)

power_8_A <- get_power(n_sim = 10000, scenario = scenario_8_A, t_star = NULL, s_star = 0.5)
power_8_B <- get_power(n_sim = 10000, scenario = scenario_8_B, t_star = NULL, s_star = 0.5)
power_8_C <- get_power(n_sim = 10000, scenario = scenario_8_C, t_star = NULL, s_star = 0.5)

power_9_A <- get_power(n_sim = 10000, scenario = scenario_9_A, t_star = NULL, s_star = 0.5)
power_9_B <- get_power(n_sim = 10000, scenario = scenario_9_B, t_star = NULL, s_star = 0.5)
power_9_C <- get_power(n_sim = 10000, scenario = scenario_9_C, t_star = NULL, s_star = 0.5)


################################################################

data_1 <- rbind(power_1_A, power_1_B, power_1_C) %>% 
  cbind(treatment_effect = rep("Scenario 1 (Homogeneous HRs)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_2 <- rbind(power_2_A, power_2_B, power_2_C) %>% 
  cbind(treatment_effect = rep("Scenario 2 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_3 <- rbind(power_3_A, power_3_B, power_3_C) %>% 
  cbind(treatment_effect = rep("Scenario 3 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")


################################################################

data_4 <- rbind(power_4_A, power_4_B, power_4_C) %>% 
  cbind(treatment_effect = rep("Scenario 4 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_5 <- rbind(power_5_A, power_5_B, power_5_C) %>% 
  cbind(treatment_effect = rep("Scenario 5 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")


data_6 <- rbind(power_6_A, power_6_B, power_6_C) %>% 
  cbind(treatment_effect = rep("Scenario 6 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

###################################################################

data_7 <- rbind(power_7_A, power_7_B, power_7_C) %>% 
  cbind(treatment_effect = rep("Scenario 7 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_8 <- rbind(power_8_A, power_8_B, power_8_C) %>% 
  cbind(treatment_effect = rep("Scenario 8 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_9 <- rbind(power_9_A, power_9_B, power_9_C) %>% 
  cbind(treatment_effect = rep("Scenario 9 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

###################################################################

null_ph_nph <- rbind(cbind(rbind(data_1, data_2, data_3),
                           pattern = "Proportional Hazards"),
                     cbind(rbind(data_4, data_5, data_6),
                           pattern = "Anticipated Delayed Effects"),
                     cbind(rbind(data_7, data_8, data_9),
                           pattern = "Null Overall Treatment Effect")) 

null_ph_nph$prognostic_effect <- factor(null_ph_nph$prognostic_effect,
                                        levels = c("None", "Moderate", "Strong"))

null_ph_nph$pattern <- factor(null_ph_nph$pattern,
                              levels = c("Proportional Hazards", "Anticipated Delayed Effects", "Null Overall Treatment Effect"))

save(null_ph_nph, file = "null_ph_nph_s_star.RData")

#Plos with all tests
plot_power = ggplot(data = null_ph_nph,
                    mapping = aes(x = prognostic_effect,
                                  y = as.numeric(as.character(value)),
                                  fill = test)) +
  facet_wrap(pattern ~ treatment_effect, scales = "free_y") +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_brewer(palette="Accent",
                    labels = c(expression("Strat. LR " (tilde(Z))), 
                               expression("Strat. LR from Mehrotra et al. " (Z^n)),
                               expression("Strat. WLR on sample size scale " (tilde(Z)^wn)),
                               expression("Strat. WLR on U-statistic scale " (tilde(Z)^wu)),
                               expression("Strat. WLR on Z-statistic scale " (tilde(Z)^wz)),
                               expression("Unstrat. LR " (Z)),
                               expression("Unstrat. WLR " (Z^w)))) +
  theme_bw() +
  #geom_hline(yintercept = 0.9, linetype = 2) + 
  geom_hline(yintercept = 0.025, linetype = 2) +
  scale_y_continuous("Power") +
  xlab("Prognostic effect") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text.align = 0)


plot_power
# Labels of the Z statistics with equations
# scale_fill_discrete(labels = c(expression(tilde(Z)^LR), 
#                                expression("Z"^n),
#                                expression(tilde(Z)^wn),
#                                expression(tilde(Z)^wu),
#                                expression(tilde(Z)^wz),
#                                expression(Z^LR),
#                                expression("Z"^w)))


ggsave("plot_tests_2_s_star.pdf", plot_power, width = 10, height = 7, dpi = 300)







