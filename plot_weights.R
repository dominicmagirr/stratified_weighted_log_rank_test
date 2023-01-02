library(tidyverse)
library(survival)
library(ggpubr)
library(survminer)
library(nphRCT)
##############################################
#OAK TRIAL
## this data can be found at https://doi.org/10.1038/s41591-018-0134-3  
data_oak <- readxl::read_excel("41591_2018_134_MOESM3_ESM.xlsx",sheet = 3) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(OS.EVENT = -1 * (OS.CNSR - 1),
         ARM = ifelse(TRT01P == "MPDL3280A","Atezolizumab","Docetaxel"))

data_oak_1 <- data_oak %>% filter(ECOGGR == 1)
data_oak_0 <- data_oak %>% filter(ECOGGR == 0)

head(data_oak)

at_risk_1 <- nphRCT::find_at_risk(Surv(OS, OS.EVENT) ~ TRT01P,
                                data = data_oak_1,
                                include_cens = FALSE)


at_risk_0 <- nphRCT::find_at_risk(Surv(OS, OS.EVENT) ~ TRT01P,
                                  data = data_oak_0,
                                  include_cens = FALSE)

w_t_1 <- nphRCT::find_weights(Surv(OS, OS.EVENT) ~ TRT01P,
                               data = data_oak_1,
                               method = "mw",
                               t_star = 12)

w_t_0 <- nphRCT::find_weights(Surv(OS, OS.EVENT) ~ TRT01P,
                                 data = data_oak_0,
                                 method = "mw",
                                 t_star = 12)


w_s_1 <- nphRCT::find_weights(Surv(OS, OS.EVENT) ~ TRT01P,
                              data = data_oak_1,
                              method = "mw",
                              s_star = 0.5)

w_s_0 <- nphRCT::find_weights(Surv(OS, OS.EVENT) ~ TRT01P,
                              data = data_oak_0,
                              method = "mw",
                              s_star = 0.5)




####################################
## find V and V_W and update weights
v_1 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_1,
                      method = "mw",
                      t_star = 0)$v_u

v_0 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_0,
                      method = "mw",
                      t_star = 0)$v_u


v_t_1 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_1,
                      method = "mw",
                      t_star = 12)$v_u

v_t_0 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_0,
                      method = "mw",
                      t_star = 12)$v_u

v_s_1 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_1,
                      method = "mw",
                      s_star = 0.5)$v_u

v_s_0 <- nphRCT::wlrt(Surv(OS, OS.EVENT) ~ TRT01P,
                      data = data_oak_0,
                      method = "mw",
                      s_star = 0.5)$v_u



###################################
## plot U,Z,N weights (t_star)
par(mfrow = c(2,3))
plot(at_risk_1$t_j, w_t_1, ylim = c(0,3),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on U-scale (t*)", type = "l")


points(at_risk_0$t_j, w_t_0, type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))

plot(at_risk_1$t_j, w_t_1 * sqrt(v_1 / v_t_1), ylim = c(0,1.5),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on Z-scale (t*)", type = "l")
points(at_risk_0$t_j, w_t_0 * sqrt(v_0 / v_t_0), type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))

plot(at_risk_1$t_j, w_t_1 * (535 / v_1), ylim = c(0,15),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on N-scale (t*)", type = "l")
points(at_risk_0$t_j, w_t_0 * (315 / v_0) , type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))


###################################
## plot U,Z,N weights (s_star)
plot(at_risk_1$t_j, w_s_1, ylim = c(0,3),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on U-scale (s*)", type = "l")


points(at_risk_0$t_j, w_s_0, type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))

plot(at_risk_1$t_j, w_s_1 * sqrt(v_1 / v_s_1), ylim = c(0,1.5),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on Z-scale (s*)", type = "l")
points(at_risk_0$t_j, w_s_0 * sqrt(v_0 / v_s_0), type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))

plot(at_risk_1$t_j, w_s_1 * (535 / v_1), ylim = c(0,15),
     xlab = "Event time",
     ylab = "Weight",
     main = "Combine on N-scale (s*)", type = "l")
points(at_risk_0$t_j, w_s_0 * (315 / v_0) , type = "l", lty = 2)
legend("bottomright",  c("ECOG = 1", "ECOG = 0"), lty = c(1,2))


