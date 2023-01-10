library(tidyverse)
library(survival)
library(ggpubr)
library(survminer)
library(nphRCT)

#This file produces Figures 1 and 2 as well as Table 2.

#######################################
#Set the output and data folders paths#
#######################################

data_path = "./Data/"
output_path = "./Outputs/"

#########################
#OAK and POPLAR trials
#########################

poplar_plot <- list()
oak_plot <- list()

#POPLAR TRIAL

data_poplar <- readxl::read_excel(paste0(data_path,"41591_2018_134_MOESM3_ESM.xlsx"),sheet = 2) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(OS.EVENT = -1 * (OS.CNSR - 1),
         ARM = ifelse(TRT01P == "MPDL3280A","Atezolizumab","Docetaxel"))

km_poplar <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_poplar)
km_poplar_strata1 <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_poplar %>% filter(ECOGGR == 1))
km_poplar_strata2 <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_poplar %>% filter(ECOGGR == 0))

poplar_plot[[1]] = survminer::ggsurvplot(km_poplar_strata1, 
                                         data = data_poplar %>% filter(ECOGGR == 1), 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "First stratum (ECOG = 1)",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = "Survival probability",
                                         risk.table.fontsize = 2.5,
                                         legend = "none")

poplar_plot[[2]] = survminer::ggsurvplot(km_poplar_strata2, 
                                         data = data_poplar %>% filter(ECOGGR == 0), 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "Second stratum (ECOG = 0)",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = " ",
                                         risk.table.fontsize = 2.5,
                                         legend = "none")

poplar_plot[[3]] = survminer::ggsurvplot(km_poplar, 
                                         data = data_poplar, 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "Overall",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = " ",
                                         risk.table.fontsize = 2.5,
                                         legend = c(0.8,0.9))

poplar_plots = arrange_ggsurvplots(poplar_plot, ncol = 3, nrow = 1)

ggsave(paste0(output_path,"Figure_2.pdf"), poplar_plots, width = 13, height = 5, dpi = 300)


#OAK TRIAL

data_oak <- readxl::read_excel(paste0(data_path,"41591_2018_134_MOESM3_ESM.xlsx"),sheet = 3) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(OS.EVENT = -1 * (OS.CNSR - 1),
         ARM = ifelse(TRT01P == "MPDL3280A","Atezolizumab","Docetaxel"))



km_oak <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_oak)
km_oak_strata1 <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_oak %>% filter(ECOGGR == 1))
km_oak_strata2 <- survfit(Surv(OS, OS.EVENT) ~ ARM, data = data_oak %>% filter(ECOGGR == 0))

oak_plot[[1]] = survminer::ggsurvplot(km_oak_strata1, 
                                         data = data_oak %>% filter(ECOGGR == 1), 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "First stratum (ECOG = 1)",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = "Survival probability",
                                         risk.table.fontsize = 2.5,
                                         legend = "none")

oak_plot[[2]] = survminer::ggsurvplot(km_oak_strata2, 
                                         data = data_oak %>% filter(ECOGGR == 0), 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "Second stratum (ECOG = 0)",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = " ",
                                         risk.table.fontsize = 2.5,
                                         legend = "none")

oak_plot[[3]] = survminer::ggsurvplot(km_oak, 
                                         data = data_oak, 
                                         legend.labs = c("Docetaxel","Atezolizumab"),
                                         risk.table = T, 
                                         conf.int = FALSE,
                                         title = "Overall",
                                         surv.median.line = "hv",
                                         break.x.by = 2,
                                         legend.title = "",
                                         xlab = "Time (months)",
                                         ylab = " ",
                                         risk.table.fontsize = 2.5,
                                         legend = c(0.8,0.9))


oak_plots = arrange_ggsurvplots(oak_plot, ncol = 3, nrow = 1)

ggsave(paste0(output_path,"Figure_1.pdf"), oak_plots, width = 13, height = 5, dpi = 300)





##########################################
#  ------   Tests from TABLE 2   ------  #
##########################################

get_table_2 <- function(dat, t_star){
  
  dat_final = dat %>%
    mutate(time = OS,
           group = ifelse(ARM == "Docetaxel", "experimental","control"),
           event = ifelse(OS.EVENT == 1, TRUE, FALSE),
           strata = ifelse(ECOGGR == 1, "first","second"))
  
  dat_final$group <- factor(dat_final$group, levels = c("experimental", "control"))
  
  
  res_strata <- wlrt(formula=Surv(time,event)~group+strata(strata),
                     data=dat_final,
                     method="mw",
                     t_star = t_star
  )
  
  
  
  res_strata_0 <- wlrt(formula=Surv(time,event)~group+strata(strata),
                       data=dat_final,
                       method="mw",
                       t_star = 0
  )
  
  res_unstrata <- wlrt(formula=Surv(time,event)~group,
                       data=dat_final,
                       method="mw",
                       t_star = t_star
  )
  
  res_unstrata_0 <- wlrt(formula=Surv(time,event)~group,
                         data=dat_final,
                         method="mw",
                         t_star = 0
  )
  
  z <- res_unstrata_0$z
  z_w <- res_unstrata$z
  tilde_z <- res_strata_0$combined$z
  tilde_z_w_u <- sum(res_strata$by_strata$u) / sqrt(sum(res_strata$by_strata$v_u))
  
  #### create tilde_z_w_z
  v_u <- res_strata_0$by_strata$v_u
  tilde_z_w_z <- sum(sqrt(v_u) * res_strata$by_strata$z) / sqrt(sum(v_u))
  
  #### create tilde_z_w_n
  v_u_w <- res_strata$by_strata$v_u
  
  n_first <- dplyr::pull(dat_final %>% filter(strata == "first") %>% summarise(n = n()))
  n_second <- dplyr::pull(dat_final %>% filter(strata == "second") %>% summarise(n = n()))
  
  u_w_n <- sum(res_strata$by_strata$u / v_u_w * c(n_first, n_second))
  v_w_n <- sum(c(n_first, n_second) ^ 2 / v_u_w)
  tilde_z_w_n <- u_w_n / sqrt(v_w_n)
  
  
  #### create tilde_z_n
  u_n <- sum(res_strata_0$by_strata$u / v_u * c(n_first, n_second))
  v_n <- sum(c(n_first, n_second) ^ 2 / v_u)
  tilde_z_n <- u_n / sqrt(v_n)
  
  
  #### return z's (these values are presented in Table 2)
  data.frame(z = z,
             z_w = z_w,
             tilde_z = tilde_z,
             tilde_z_n = tilde_z_n,
             tilde_z_w_u = tilde_z_w_u,
             tilde_z_w_z = tilde_z_w_z,
             tilde_z_w_n = tilde_z_w_n) %>% round(2) 
  
}

table_2 <- rbind(get_table_2(dat = data_oak, t_star = 6),
                 get_table_2(dat = data_oak, t_star = 12),
                 get_table_2(dat = data_poplar, t_star = 6),
                 get_table_2(dat = data_poplar, t_star = 12))

rownames(table_2) = c("OAK, t* = 6",
                     "OAK, t* = 12",
                     "POPLAR, t* = 6",
                     "POPLAR, t* = 12")



save(table_2, file = paste0(output_path, "table_2.RData"))






