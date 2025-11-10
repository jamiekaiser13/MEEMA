library(tidyverse)
library(ggplot2)

#BMI
plot_BMI <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=BMI, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  ylab("Body Mass Index") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(BMI, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(limits= c(20, 55), expand = c(0.15, 0))

BMI_change_intervention <- c(0.2, -1.5, -2.5, -1.8)
mean(BMI_change_intervention)
#[1] -1.4
BMI_change_control <- c(0.4)
t.test(BMI_change_control, BMI_change_intervention, alternative = "greater", var.equal = TRUE)
#Not sig, but we shouldn't use a t-test anyway with this small sample size because we can't accurately calculate variance for control, not normally distributed.

wilcox.test(BMI_change_control, BMI_change_intervention)
# Wilcoxon rank sum exact test
# 
# data:  BMI_change_control and BMI_change_intervention
# W = 4, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0

#trying to use a one-sample Wilcox test against a hypothesized median, in which
#the median is the change in BMI for the control individual (0.4)
#alternative hypothesis: the change in BMI for the intervention group is less than
#the change for the individual in the control group
wilcox.test(BMI_change_intervention, mu = 0.4, alternative = "less")

#Percent fat mass
plot_percentfat <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=PercentFat, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  ylab("Percent Fat Mass") +
  theme_minimal() + 
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(PercentFat, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

#Lean mass
plot_leanmass <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=LM, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  ylab("Lean Mass in g") +
  theme_minimal() +
  guides(linetype = guide_legend(title = "Group")) +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(LM, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

meema_metadata_moms$SampleID
meema_metadata_moms$LM
leanmass <- data.frame(SampleID= meema_metadata_moms$SampleID, 
                       LM= meema_metadata_moms$LM)
leanmass
34498-34290
48232-48961
46367-44949
85501-35489
45170-40343
lm_delta_ctrl <- 4827
lm_delat_intv <- c(208, -729, 1418, 50012)
wilcox.test(lm_delta_ctrl, lm_delat_intv)
# Wilcoxon rank sum exact test
# 
# data:  lm_delta_ctrl and lm_delat_intv
# W = 3, p-value = 0.8
# alternative hypothesis: true location shift is not equal to 0


#Bone mineral content
plot_BMC <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=BMC, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(linetype = guide_legend(title = "Group")) +
  ylab("Bone Mineral Content in g") + 
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(BMC, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

bmc <- data.frame(SampleID= meema_metadata_moms$SampleID, 
                       BMC= meema_metadata_moms$BMC)
bmc
2107-2181
2779-2811
2829-2893
2061-2127
3122-3415
bmc_delta_ctrl <- -293
bmc_delta_intv <- c(-74, -32, -64, -66)
wilcox.test(bmc_delta_ctrl, bmc_delta_intv)
# Wilcoxon rank sum exact test
# 
# data:  bmc_delta_ctrl and bmc_delta_intv
# W = 0, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0
median(bmc_delta_intv)


#Bone mineral density
plot_BMD <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=BMD, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(linetype = guide_legend(title = "Group")) +
  ylab(expression(paste("Bone Mineral Density in g/cm"^"2"))) + 
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(BMD, 3)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

bmd <- data.frame(SampleID= meema_metadata_moms$SampleID, 
                  BMD= meema_metadata_moms$BMD)
bmd
1.083-1.092
1.215-1.230
1.256-1.260
1.088-1.102
1.317-1.341
bmd_delta_ctrl <- -0.024
bmd_delta_intv <- c(-0.009, -0.015, -0.004, -0.014)
wilcox.test(bmd_delta_ctrl, bmd_delta_intv)
# Wilcoxon rank sum exact test
# 
# data:  bmd_delta_ctrl and bmd_delta_intv
# W = 0, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0
median(bmd_delta_intv)

#HEI score
plot_HEI <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=HEI2015.TOTAL.SCORE, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(linetype = guide_legend(title = "Group")) +
  ylab("HEI score") + 
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(HEI2015.TOTAL.SCORE, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))


#Fat Free Mass
plot_fatfree <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=FFM, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(linetype = guide_legend(title = "Group")) +
  ylab("Fat Free Mass in g") + 
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(FFM, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))


#VO2max
plot_VO2max <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=VO2max, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2, show.legend = FALSE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Predicted VO2 Max (mL/kg/min)") + 
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(VO2max, 1)), 
            hjust= (ifelse(meema_metadata_moms$Timepoint == "Baseline", 1.5, -0.5)), size=2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

#plot for the manuscript
plot_VO2max_man <- ggplot(meema_metadata_moms, aes(x=Timepoint, y=VO2max, color=Group, group = Baseline.cpartid)) + 
  geom_line(aes(linetype=Group), size = 1.2) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Predicted VO2 Max (mL/kg/min)") + 
  xlab("") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  scale_x_discrete(expand = c(0.08, 0)) + scale_y_continuous(expand = c(0.15, 0))
plot_VO2max_man

VO2Max <- data.frame(SampleID= meema_metadata_moms$SampleID, 
                  VO2max= meema_metadata_moms$VO2max)
VO2Max
25.24-23.33
30.24-28.56
30.59-30.28
28.85-24.64
24.52-26.22
VO2Max_delta_ctrl <- -1.7
VO2Max_delta_intv <- c(1.91, 1.68, 0.31, 4.21)
wilcox.test(VO2Max_delta_ctrl, VO2Max_delta_intv)
# Wilcoxon rank sum exact test
# 
# data:  VO2Max_delta_ctrl and VO2Max_delta_intv
# W = 0, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0

#wilcoxon signed rank 
wilcox.test(VO2Max_delta_intv, mu = -1.7, alternative = "greater")


mean(VO2Max_delta_intv)
# [1] 2.0275
median(VO2Max_delta_intv)
# [1] 1.795

#Vitamin D
#there is dietary metadata that I'm going to pull from that gives the vitamin D 
#consumption
vitaminD_calcium_metadata <- read.csv("C:/Users/jmktw/Downloads/MEEMA-VitaminD-Calcium.csv", stringsAsFactors = TRUE)
#DSVDBA
DSVDB <- ggplot(vitaminD_calcium_metadata, aes(x=Timepoint, y=DSVDBA, color= Group, group = cpartid)) +
  geom_line(aes(linetype = Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Vitamin D Intake from Dietary Supplements in mcg") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(DSVDBA, 1)),
            hjust = (ifelse(vitaminD_calcium_metadata$Timepoint == "Baseline", 1.5, -0.5)),
            size = 2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))
#NTVDBA
NTVDBA <- ggplot(vitaminD_calcium_metadata, aes(x=Timepoint, y=NTVDBA, color= Group, group = cpartid)) +
  geom_line(aes(linetype = Group), linewidth = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("NTVDBA") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(NTVDBA, 1)),
            hjust = (ifelse(vitaminD_calcium_metadata$Timepoint == "Baseline", 1.5, -0.5)),
            size = 2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))
#RDVitDbA
plot_RDVitDbA <- ggplot(vitaminD_calcium_metadata, aes(x=Timepoint, y=RDVitDbA, color= Group, group = cpartid)) +
  geom_line(aes(linetype = Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Percent Recommended Daily Value of Vitamin D") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(RDVitDbA, 1)),
            hjust = (ifelse(vitaminD_calcium_metadata$Timepoint == "Baseline", 1.5, -0.5)),
            size = 2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))
#DSCBA
plot_DSCBA <- ggplot(vitaminD_calcium_metadata, aes(x=Timepoint, y=DSCBA, color= Group, group = cpartid)) +
  geom_line(aes(linetype = Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Calcium Intake from Dietary Supplements in mg") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(DSCBA, 1)),
            hjust = (ifelse(vitaminD_calcium_metadata$Timepoint == "Baseline", 1.5, -0.5)),
            size = 2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))
#RDCalbA
plot_RDCalbA <- ggplot(vitaminD_calcium_metadata, aes(x=Timepoint, y=RDCalbA, color= Group, group = cpartid)) +
  geom_line(aes(linetype = Group), size = 1.2, show.legend = TRUE) +
  geom_point(size = 3.5) +
  guides(color = guide_legend(title = "Group")) +
  ylab("Percent Recommended Daily Value of Calcium") +
  theme_minimal() +
  scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  geom_text(aes(label=round(RDCalbA, 1)),
            hjust = (ifelse(vitaminD_calcium_metadata$Timepoint == "Baseline", 1.5, -0.5)),
            size = 2.6) +
  scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0))

onedrive <- "C:/Users/jmktw/OneDrive - University of Illinois Chicago/MEEMA Paper"
setwd(onedrive)
delta_diversity_metrics <- read.csv("Delta_Metrics.csv")

View(delta_metrics)



ctrl_delta_fat <- -4.2
intv_delta_fat <- c(0.2, -3.9, -4, -0.7)
wilcox.test(ctrl_delta_fat, intv_delta_fat)
# Wilcoxon rank sum exact test
# 
# data:  ctrl_delta_fat and intv_delta_fat
# W = 0, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0


delta_VO2max_df
intv_VO2_delta <- c(1.91, 1.68, 0.31, 4.21)
ctrl_VO2_delta <- -1.7

wilcox.test(ctrl_VO2_delta,intv_VO2_delta)
# Wilcoxon rank sum exact test
# 
# data:  ctrl_VO2_delta and intv_VO2_delta
# W = 4, p-value = 0.4
# alternative hypothesis: true location shift is not equal to 0

