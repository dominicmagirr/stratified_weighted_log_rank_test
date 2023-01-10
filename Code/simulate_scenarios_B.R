library(tidyr)
library(nphRCT)
#############################################
#Set the output, data and code folders paths#
#############################################

data_path = "./Data/"
output_path = "./Outputs/"

#Produce Figure 4
load(paste0(data_path, "null_ph_nph_t_star.RData"))

plot_figure_4 = ggplot(data = null_ph_nph,
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
  geom_hline(yintercept = 0.025, linetype = 2) +
  scale_y_continuous("Power") +
  xlab("Prognostic effect") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text.align = 0)


plot_figure_4


ggsave(paste0(output_path,"Figure_4.pdf"), plot_figure_4, width = 10, height = 7, dpi = 300)


#Produce Figure 5
load(paste0(data_path, "null_ph_nph_s_star.RData"))

plot_figure_5 = ggplot(data = null_ph_nph,
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
  geom_hline(yintercept = 0.025, linetype = 2) +
  scale_y_continuous("Power") +
  xlab("Prognostic effect") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text.align = 0)


plot_figure_5

ggsave(paste0(output_path,"Figure_5.pdf"), plot_figure_5, width = 10, height = 7, dpi = 300)

