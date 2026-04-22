library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
library(survminer)
library(grid)
thedata=read.csv(here("Fitting of Natural History Model/Jane/toy_analysis.csv"))%>%
  mutate(combo=paste0(OMST,"_",LMST,"_",Sens.early....,"_",Sens.late....,"_",Specificity....),
         design=paste0(X..screens,"_",Follow.up ))


n_plot_traditional = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90 & LMST==1),
  aes(x = design, y = N.Traditional)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  ylim(0, 125000) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("N Traditional") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")

n_plot_IE = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90&LMST==1),
  aes(x = design, y = N.IE)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  ylim(0, 125000) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("N IE") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")

n_plot_PA = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90&LMST==1),
  aes(x = design, y = N.PA)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  ylim(0, 125000) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("N PA") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")

RE_plot_PA = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90&LMST==1),
  aes(x = design, y = N_Traditional...N_PA)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("RE PA") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")

RE_plot_IE = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90&LMST==1),
  aes(x = design, y = N_Traditional...N_IE)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("RE IE") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")

sample_size_RR = ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90 &Follow.up==2 & X..screens==3),
  aes(x = RR.Traditional , y = N.IE)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ggtitle("Sample size versus effect size") +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity")





################################################
#Chat GPT combined plot
###############################################
library(ggplot2)
library(cowplot)
library(grid)

# Remove legends from all plots
p1 <- n_plot_traditional + theme(legend.position = "none")
p2 <- n_plot_IE + theme(legend.position = "none")
p3 <- n_plot_PA + theme(legend.position = "none")
p4 <- RE_plot_PA + theme(legend.position = "none")
p5 <- RE_plot_IE + theme(legend.position = "none")

legend_plot <- ggplot(
  data = subset(thedata, Specificity.... == 100 & Sens.late.... == 90&LMST==1),
  aes(x = design, y = N.Traditional)
) +
  geom_point(aes(color = combo), shape = 1) +
  geom_line(aes(color = combo, group = combo)) +
  theme_bw() +
  labs(colour = "OMST_LMST_EarlySens_LateSens_Specificity") +
  theme(legend.position = "right")

legend <- get_legend(legend_plot)

legend_panel <- ggdraw() +
  draw_grob(legend, x = 0.05, y = 0.05, width = 0.9, height = 0.9)

panel_grid <- plot_grid(
  p1, p2, p3,
  legend_panel, p5, p4,
  ncol = 3,
  align = "hv"
)

xlab_grob <- textGrob(
  "No. screens _ Years of follow-up after last screen",
  gp = gpar(fontsize = 11)
)

final_plot <- plot_grid(
  panel_grid,
  ggdraw() + draw_grob(xlab_grob),
  ncol = 1,
  rel_heights = c(1, 0.06)
)

final_plot