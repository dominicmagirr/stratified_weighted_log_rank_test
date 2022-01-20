

#########################################################
#########################################################
## Recruitment
#########################################################
#########################################################
recruitment_first  = list(n_0 = 86, n_1 = 86, r_period = 9, k = 1)
recruitment_second = list(n_0 = 86, n_1 = 86, r_period = 9, k = 1)

##########################################################
##########################################################
## Models
##########################################################
##########################################################

##########################################
## 0. Null effect
##########################################

###### A. Non-prognostic #################

model_0_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 1))

model_0_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 1))

###### B. moderate prognostic ############


model_0_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 1, log(2) / 6 * 1))

model_0_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 1, log(2) / 10 * 1))

###### C. strong prognostic ##############


model_0_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 1, log(2) / 3 * 1))

model_0_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 1, log(2) / 15 * 1))

##########################################
## 1. Homogeneous
##########################################

###### A. Non-prognostic #################

model_1_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 2/3, log(2) / 8 * 2/3))

model_1_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 2/3, log(2) / 8 * 2/3))

###### B. moderate prognostic ############


model_1_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 2/3, log(2) / 6 * 2/3))

model_1_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 2/3, log(2) / 10 * 2/3))

###### C. strong prognostic ##############


model_1_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 2/3, log(2) / 3 * 2/3))

model_1_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 2/3, log(2) / 15 * 2/3))


##########################################
## 2. Quantitative difference (poor better)
##########################################

###### A. Non-prognostic #################

model_2_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 0.6, log(2) / 8 * 0.6))

model_2_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 3/4, log(2) / 8 * 3/4))

###### B. moderate prognostic ############


model_2_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 0.6, log(2) / 6 * 0.6))

model_2_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 3/4, log(2) / 10 * 3/4))

###### C. strong prognostic ##############


model_2_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 0.6, log(2) / 3 * 0.6))

model_2_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 3/4, log(2) / 15 * 3/4))


##########################################
## 3. Quantitative difference (poor worse)
##########################################

###### A. Non-prognostic #################

model_3_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 3/4, log(2) / 8 * 3/4))

model_3_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 0.6, log(2) / 8 * 0.6))

###### B. moderate prognostic ############


model_3_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 3/4, log(2) / 6 * 3/4))

model_3_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 0.6, log(2) / 10 * 0.6))

###### C. strong prognostic ##############


model_3_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 3/4, log(2) / 3 * 3/4))

model_3_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 0.6, log(2) / 15 * 0.6))

##########################################
## 4. Qualitative difference (poor better)
##########################################

###### A. Non-prognostic #################

model_4_A_first = list(change_points = c(11),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 2/3, log(2) / 8 * 0.72))

model_4_A_second = list(change_points = c(11),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 1.5, log(2) / 8 * 2.2))

###### B. moderate prognostic ############


model_4_B_first = list(change_points = c(11),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 2/3, log(2) / 6 * 2/3))

model_4_B_second = list(change_points = c(11),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 1.5, log(2) / 10 * 1.5))

###### C. strong prognostic ##############


model_4_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 2/3, log(2) / 3 * 0.3))

model_4_C_second = list(change_points = c(8),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 2.1, log(2) / 15 * 1.2))

##########################################
## 5. Qualitative difference (poor worse)
##########################################

###### A. Non-prognostic #################

model_5_A_first = list(change_points = c(9),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 1.8, log(2) / 8 * 2))

model_5_A_second = list(change_points = c(9),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 0.5, log(2) / 8 * 0.8))

###### B. moderate prognostic ############


model_5_B_first = list(change_points = c(11),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 1.5, log(2) / 6 * 2.4))

model_5_B_second = list(change_points = c(11),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 0.65, log(2) / 10 * 0.95))

###### C. strong prognostic ##############


model_5_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 1.4, log(2) / 3 * 2))

model_5_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 0.675, log(2) / 15 * 1.075))

##########################################
## 6. "Homogeneous"
##########################################

###### A. Non-prognostic #################

model_6_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 1/2))

model_6_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 1/2))

###### B. moderate prognostic ############


model_6_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 1, log(2) / 6 * 0.4))

model_6_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 1, log(2) / 10 * 0.5))

###### C. strong prognostic ##############


model_6_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 1, log(2) / 3 * 0.2))

model_6_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 1, log(2) / 15 * 0.5))


##########################################
## 7. Quantitative difference (poor better)
##########################################

###### A. Non-prognostic #################

model_7_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 0.3))

model_7_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 0.7))

###### B. moderate prognostic ############


model_7_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 1, log(2) / 6 * 0.3))

model_7_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 1, log(2) / 10 * 0.7))

###### C. strong prognostic ##############


model_7_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 1, log(2) / 3 * 0.1))

model_7_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 1, log(2) / 15 * 0.8))


##########################################
## 8. Quantitative difference (poor worse)
##########################################

###### A. Non-prognostic #################

model_8_A_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 8, log(2) / 8),
                       lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 0.7))

model_8_A_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 8, log(2) / 8),
                        lambdas_1 = c(log(2) / 8 * 1, log(2) / 8 * 0.3))

###### B. moderate prognostic ############


model_8_B_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 6, log(2) / 6),
                       lambdas_1 = c(log(2) / 6 * 1, log(2) / 6 * 0.7))

model_8_B_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 10, log(2) / 10),
                        lambdas_1 = c(log(2) / 10 * 1, log(2) / 10 * 0.3))

###### C. strong prognostic ##############


model_8_C_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 3, log(2) / 3),
                       lambdas_1 = c(log(2) / 3 * 1, log(2) / 3 * 0.7))

model_8_C_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 1, log(2) / 15 * 0.3))


##########################################
## OAK
##########################################


model_OAK_first = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 7.5, log(2) / 7.5),
                       lambdas_1 = c(log(2) / 7.5, log(2) / 13))

model_OAK_second = list(change_points = c(6),
                        lambdas_0 = c(log(2) / 15, log(2) / 15),
                        lambdas_1 = c(log(2) / 15 * 1, log(2) / 17.5))
