library(tidyverse)
library(survival)
library(ggpubr)
#devtools::load_all("~/NPH/swlrt") # not needed now using swlrt_fast.R

source("set_up_model.R")
source("swlrt_fast.R")
source("wlrt_fast.R")

####################################################
####################################################
## simulation (as a function)
####################################################
####################################################

###########################################################
####### simulate from a piece-wise exponential distribution
t_piecewise_exp <- function(n = 10, 
                            change_points = c(6, 12),
                            lambdas = c(log(2) / 9, log(2) / 9, log(2) / 9)){
  
  t_lim <- matrix(rep(c(diff(c(0, change_points)), Inf), each = n), nrow = n)
  t_sep <- do.call(cbind, purrr::map(lambdas, rexp, n = n))
  which_cells <- t(apply(t_sep < t_lim, 1, function(x){
    rep(c(T,F), c(min(which(x)), length(x) - min(which(x))))
  } ))
  rowSums(pmin(t_sep, t_lim) * which_cells)
}
##########################################################
######## simulate data according to trial parameters
sim_t_uncensored <- function(model,
                             recruitment){
  
  rec_0 <- recruitment$r_period * runif(recruitment$n_0) ^ (1 / recruitment$k)
  rec_1 <- recruitment$r_period * runif(recruitment$n_1) ^ (1 / recruitment$k)
  
  time_0 <- t_piecewise_exp(recruitment$n_0, model$change_points, model$lambdas_0)
  time_1 <- t_piecewise_exp(recruitment$n_1, model$change_points, model$lambdas_1)
  
  data.frame(time = c(time_0, time_1),
             rec = c(rec_0, rec_1),
             group = rep(c("control", "experimental"), c(recruitment$n_0, recruitment$n_1)))
  
}
#########################################################
## apply data cut off to simulated data set
apply_dco <- function(df,
                      dco = NULL,
                      events = NULL){
  
  if (is.null(dco) && is.null(events)) stop("Must specify either dco or events")
  
  df$cal_time <- df$time + df$rec
  
  if (is.null(dco)){
    dco <- sort(df$cal_time)[events]
  }
  
  df_dco <- df[df$rec < dco, ]
  df_dco$event <-  df_dco$cal_time <= dco
  df_dco$time <- pmin(df_dco$time, dco - df_dco$rec)
  df_dco$dco <- dco
  
  df_dco
  
}

########################################################
# model_first <- model_0_A_first
# model_second <- model_0_A_second
# 
# recruitment_first <- recruitment_first
# recruitment_second <- recruitment_second
# 
# dco <- 24
# 
# t_star <- 12
# 
# 
# sim_1_trial(dummy, 
#             model_first, 
#             model_second,
#             recruitment_first,
#             recruitment_second,
#             dco,
#             t_star)
########################################################
sim_1_trial <- function(dummy, 
                        model_first, 
                        model_second,
                        recruitment_first,
                        recruitment_second,
                        dco,
                        t_star){
  
  df_uncensored_first <- sim_t_uncensored(model_first, recruitment_first)
  df_uncensored_second <- sim_t_uncensored(model_second, recruitment_second)
  
  df_final_first <- apply_dco(df_uncensored_first, dco = dco)
  df_final_second <- apply_dco(df_uncensored_second, dco = dco)
  
  df_final_first$strata <- "first"
  df_final_second$strata <- "second"
  
  df_final <- rbind(df_final_first, df_final_second)
  
  #########################################
  ### apply swlrt
  res_strata <- swlrt_fast(df = df_final,
                           trt_colname = "group",
                           time_colname = "time",
                           event_colname = "event",
                           strat_colname = "strata",
                           wlr = "mw",
                           t_star = t_star)
  
  res_strata_0 <- swlrt_fast(df = df_final,
                             trt_colname = "group",
                             time_colname = "time",
                             event_colname = "event",
                             strat_colname = "strata",
                             wlr = "mw",
                             t_star = 0)
  
  res_unstrata <- wlrt_fast(df = df_final,
                            trt_colname = "group",
                            time_colname = "time",
                            event_colname = "event",
                            wlr = "mw",
                            t_star = t_star)
  
  res_unstrata_0 <- wlrt_fast(df = df_final,
                              trt_colname = "group",
                              time_colname = "time",
                              event_colname = "event",
                              wlr = "mw",
                              t_star = 0)
  
  z <- res_unstrata_0$z
  z_w <- res_unstrata$z
  tilde_z <- res_strata_0$z
  tilde_z_w_u <- res_strata$z
  
  #### create tilde_z_w_z
  v_u <- res_strata_0$by_strata$v_u
  tilde_z_w_z <- sum(sqrt(v_u) * res_strata$by_strata$z) / sqrt(sum(v_u))
  
  #### create tilde_z_w_n
  v_u_w <- res_strata$by_strata$v_u
  
  n_first <- sum(df_final_first$rec < df_final_first$dco)
  n_second <- sum(df_final_second$rec < df_final_second$dco)
  
  u_w_n <- sum(res_strata$by_strata$u / v_u_w * c(n_first, n_second))
  v_w_n <- sum(c(n_first, n_second) ^ 2 / v_u_w)
  tilde_z_w_n <- u_w_n / sqrt(v_w_n)
  
  
  #### create tilde_z_n
  u_n <- sum(res_strata_0$by_strata$u / v_u * c(n_first, n_second))
  v_n <- sum(c(n_first, n_second) ^ 2 / v_u)
  tilde_z_n <- u_n / sqrt(v_n)
  
  
  #### return z's
  data.frame(z = z,
             z_w = z_w,
             tilde_z = tilde_z,
             tilde_z_w_u = tilde_z_w_u,
             tilde_z_w_z = tilde_z_w_z,
             tilde_z_n = tilde_z_n,
             tilde_z_w_n = tilde_z_w_n)
  
  
}




###########################################################################
get_power <- function(model_first, 
                      model_second,
                      recruitment_first,
                      recruitment_second,
                      dco,
                      t_star){
  
  sim_res <- purrr::map_df(1:10000, sim_1_trial, model_first, model_second, recruitment_first, recruitment_second, dco, t_star)
  apply(sim_res, 2, function(x) mean(x < qnorm(0.025)))
  
  
}
###################
set.seed(124543)
###################
power_0_A <- get_power(model_0_A_first, model_0_A_second, recruitment_first, recruitment_second, 24, 12)
power_0_B <- get_power(model_0_B_first, model_0_B_second, recruitment_first, recruitment_second, 24, 12)
power_0_C <- get_power(model_0_C_first, model_0_C_second, recruitment_first, recruitment_second, 24, 12)

power_1_A <- get_power(model_1_A_first, model_1_A_second, recruitment_first, recruitment_second, 24, 12)
power_1_B <- get_power(model_1_B_first, model_1_B_second, recruitment_first, recruitment_second, 24, 12)
power_1_C <- get_power(model_1_C_first, model_1_C_second, recruitment_first, recruitment_second, 24, 12)

power_2_A <- get_power(model_2_A_first, model_2_A_second, recruitment_first, recruitment_second, 24, 12)
power_2_B <- get_power(model_2_B_first, model_2_B_second, recruitment_first, recruitment_second, 24, 12)
power_2_C <- get_power(model_2_C_first, model_2_C_second, recruitment_first, recruitment_second, 24, 12)

power_3_A <- get_power(model_3_A_first, model_3_A_second, recruitment_first, recruitment_second, 24, 12)
power_3_B <- get_power(model_3_B_first, model_3_B_second, recruitment_first, recruitment_second, 24, 12)
power_3_C <- get_power(model_3_C_first, model_3_C_second, recruitment_first, recruitment_second, 24, 12)

power_4_A <- get_power(model_4_A_first, model_4_A_second, recruitment_first, recruitment_second, 24, 12)
power_4_B <- get_power(model_4_B_first, model_4_B_second, recruitment_first, recruitment_second, 24, 12)
power_4_C <- get_power(model_4_C_first, model_4_C_second, recruitment_first, recruitment_second, 24, 12)

power_5_A <- get_power(model_5_A_first, model_5_A_second, recruitment_first, recruitment_second, 24, 12)
power_5_B <- get_power(model_5_B_first, model_5_B_second, recruitment_first, recruitment_second, 24, 12)
power_5_C <- get_power(model_5_C_first, model_5_C_second, recruitment_first, recruitment_second, 24, 12)

power_6_A <- get_power(model_6_A_first, model_6_A_second, recruitment_first, recruitment_second, 24, 12)
power_6_B <- get_power(model_6_B_first, model_6_B_second, recruitment_first, recruitment_second, 24, 12)
power_6_C <- get_power(model_6_C_first, model_6_C_second, recruitment_first, recruitment_second, 24, 12)

power_7_A <- get_power(model_7_A_first, model_7_A_second, recruitment_first, recruitment_second, 24, 12)
power_7_B <- get_power(model_7_B_first, model_7_B_second, recruitment_first, recruitment_second, 24, 12)
power_7_C <- get_power(model_7_C_first, model_7_C_second, recruitment_first, recruitment_second, 24, 12)

power_8_A <- get_power(model_8_A_first, model_8_A_second, recruitment_first, recruitment_second, 24, 12)
power_8_B <- get_power(model_8_B_first, model_8_B_second, recruitment_first, recruitment_second, 24, 12)
power_8_C <- get_power(model_8_C_first, model_8_C_second, recruitment_first, recruitment_second, 24, 12)



################################################################
data_0 <- rbind(power_0_A, power_0_B, power_0_C) %>% 
  cbind(treatment_effect = rep("Scenario 7 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_4 <- rbind(power_4_A, power_4_B, power_4_C) %>% 
  cbind(treatment_effect = rep("Scenario 8 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_5 <- rbind(power_5_A, power_5_B, power_5_C) %>% 
  cbind(treatment_effect = rep("Scenario 9 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")
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
data_6 <- rbind(power_6_A, power_6_B, power_6_C) %>% 
  cbind(treatment_effect = rep("Scenario 4 (Homogeneous)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_7 <- rbind(power_7_A, power_7_B, power_7_C) %>% 
  cbind(treatment_effect = rep("Scenario 5 (poor prog. '>' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")

data_8 <- rbind(power_8_A, power_8_B, power_8_C) %>% 
  cbind(treatment_effect = rep("Scenario 6 (poor prog. '<' good prog.)", 3),
        prognostic_effect = c("None", "Moderate", "Strong")) %>% 
  as.data.frame() %>% 
  pivot_longer(z:tilde_z_w_n, names_to = "test")
###################################################################

null_ph_nph <- rbind(cbind(rbind(data_1, data_2, data_3),
                           pattern = "Proportional Hazards"),
                     cbind(rbind(data_6, data_7, data_8),
                           pattern = "Anticipated Delayed Effects"),
                     cbind(rbind(data_0, data_4, data_5),
                           pattern = "Null Overall Treatment Effect")) 

null_ph_nph$prognostic_effect <- factor(null_ph_nph$prognostic_effect,
                                   levels = c("None", "Moderate", "Strong"))



#Plos with all tests
plot_power = ggplot(data = null_ph_nph,
       mapping = aes(x = prognostic_effect,
                     y = as.numeric(as.character(value)),
                     fill = test)) +
  facet_wrap(pattern ~ treatment_effect) +
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
  geom_hline(yintercept = 0.9, linetype = 2) + 
  geom_hline(yintercept = 0.025, linetype = 2) +
  scale_y_continuous("Power") +
  xlab("Prognostic effect") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text.align = 0)


  # Labels of the Z statistics with equations
  # scale_fill_discrete(labels = c(expression(tilde(Z)^LR), 
  #                                expression("Z"^n),
  #                                expression(tilde(Z)^wn),
  #                                expression(tilde(Z)^wu),
  #                                expression(tilde(Z)^wz),
  #                                expression(Z^LR),
  #                                expression("Z"^w)))


ggsave("plot_tests_2.pdf", plot_power, width = 10, height = 7, dpi = 300)

###########################################################################
