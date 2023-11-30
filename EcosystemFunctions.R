### TOPTRONS ECOSYSTEM FUNCTIONS ###
## Antonia Ahme, 14.11.2023 ##

#### HOUSEKEEPING ####
require(dplyr)
require(readxl)
require(ggpubr)
require(rstatix)
require(ggplot2)
require(writexl)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/TopTrons") # set to your working directory

set.seed(22)

#### FUNCTIONS & LAYOUT ####
# Create the design for plotting
plot.theme <- theme(panel.background = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    panel.grid = element_blank(),
                    axis.title.x = element_text(size = 15, face= "bold"),
                    axis.text.x = element_text(size = 15, face= "bold"),
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_text(size = 15, face= "bold"),
                    axis.text.y = element_text(face="plain", color="black", size=15),
                    axis.ticks.y = element_blank(),
                    plot.title = element_text(size = 15, face= "bold"),
                    strip.text = element_text(size = 15, face = "bold"),
                    strip.background = element_blank(), strip.placement = "outside",
                    text = element_text(size = 15, face = "bold"), legend.position = "none")

# Color palette
temp_pal <- c("6"="skyblue2", "12"="goldenrod2", "18"="tomato2")
temp_pal2 <- c("lightskyblue1", "gold", "salmon1")

#### UPLOAD AND TIDY UP DATA ####
pan <- readxl::read_excel("Data/ParticulateNutrients.xlsx")
pan$temp <- as.factor(pan$temp)
pan_all <- pan
pan <- subset(pan, day > 12)
pan612 <- subset(pan, temp != "18")
pan618 <- subset(pan, temp != "12")

oxy <- readxl::read_excel("Data/GOP.xlsx")
oxy$temp <- as.factor(oxy$temp)
oxy_all <- oxy
oxy <- subset(oxy, day > 12)
oxy2 <- oxy
oxy$temp <- factor(oxy$temp, levels = c("6", "12", "18"))
oxy612 <- subset(oxy, temp != "18")
oxy618 <- subset(oxy, temp != "12")

#### PLOTTING: ECOSYSTEM FUNCTIONS ####
# Biomass
p1 <- ggplot(pan, aes(x=day, y=µgL_c, color=temp, group=interaction(temp,day))) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 2, position = position_dodge(width=2.25))+
  geom_point(position = position_dodge(width=2.25))+
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("POC (µg/L)")) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 4000, 500))+
  scale_color_manual(values=temp_pal)

p1

ggsave("Output/Biomass.png", p1, dpi = 300, width = 8, height = 4)

# GOP per POC
p2 <- ggplot(oxy, aes(x=day, y=o2_poc, color=temp, group=interaction(temp,day))) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 2, position = position_dodge(width=2.25))+
  geom_point(position = position_dodge(width=2.25))+
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("GOP (µmol O2 mg POC^-1 d^-1)")) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 30, 5))+
  scale_color_manual(values=temp_pal)

p2

ggsave("Output/GOP.png", p2, dpi = 300, width = 8, height = 4)

# C:N ratio
p3 <- ggplot(pan, aes(x=day, y=cn, color=temp, group=interaction(temp,day))) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 2, position = position_dodge(width=2.25))+
  geom_point(position = position_dodge(width=2.25))+
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("POC:PON (mol:mol)")) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 60, 5))+
  scale_color_manual(values=temp_pal)

p3

ggsave("Output/CN.png", p3, dpi = 300, width = 8, height = 4)

# C:P ratio
p4 <- ggplot(pan, aes(x=day, y=cp, color=temp, group=interaction(temp,day))) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 2, position = position_dodge(width=2.25))+
  geom_point(position = position_dodge(width=2.25))+
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("POC:POP (mol:mol)")) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 1100, 100))+
  scale_color_manual(values=temp_pal)

p4

ggsave("Output/CP.png", p4, dpi = 300, width = 8, height = 4)


#### STATISTICS ####
### TWO-WAY REPEATED MEASURES ANOVA
#https://stats.stackexchange.com/questions/181563/analyzing-repeated-measures-experiment-with-multiple-treatment-groups-and-multip
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
#https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

### Check bloom differences during acclimation and experiment
pan_accli <-  subset(pan_all, day < 13)

## Acclimation phase
# Check assumptions
pan2 <- pan_accli %>% select(chl_µgL, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)
pan2 <- subset(pan2, day != "0")

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(chl_µgL, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "chl_µgL",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "chl_µgL", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

XX <- subset(pan2, day == "9" & temp == "18")
ggqqplot(XX, "chl_µgL")
# looks okay

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = chl_µgL, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity not violated


## Experiment phase
# Check assumptions
pan2 <- pan %>% select(chl_µgL, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)
#pan2 <- subset(pan2, temp != "6")

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(chl_µgL, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "chl_µgL",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "chl_µgL", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = chl_µgL, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# significant temp effect and without 6 significant time effect


### Ecosystem functions
## Create a table for statistical outout of two-way RM anova
Parameter <- rep(c("Biomass", "GOP", "C:N", "C:P"), each=3)
Effect <- rep(c("Temperature","Time","Temperature:Time"), times=4)
stat <- data.frame(Parameter, Effect)
stat['DFn'] <- NA
stat['DFd'] <- NA
stat['F'] <- NA
stat['P'] <- NA

## Biomass
# Check assumptions
pan2 <- pan %>% select(µgL_c, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(µgL_c, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "µgL_c",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "µgL_c", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = µgL_c, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity not violated

res.aov_bm <- get_anova_table(res.aov)
# no significant interaction + significant effects of day and temp

# add main results to stat table
stat[1:3,3:6] <- res.aov_bm[,2:5]

# Considering no significant interaction
# comparisons for treatment variable
pan2 %>%
  pairwise_t_test(
    µgL_c ~ temp, paired = FALSE, 
    p.adjust.method = "bonferroni")

# comparisons for time variable
pan2 %>%
  pairwise_t_test(
    µgL_c ~ day, paired = TRUE, 
    p.adjust.method = "bonferroni")

## Oxygen - GOP
# Check assumptions
pan2 <- oxy2 %>% select(o2_poc, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(o2_poc, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "o2_poc",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "o2_poc", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but looks all right

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = o2_poc, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# assumption of sphericity not met - GG correction applied

res.aov_gpp <- get_anova_table(res.aov)
# no significant interaction + significant effects of day and temp

# add main results to stat table
stat[4:6,3:6] <- res.aov_gpp[,2:5]

# Considering no significant interaction
# comparisons for treatment variable
pan2 %>%
  pairwise_t_test(
    o2_poc ~ temp, paired = FALSE, 
    p.adjust.method = "bonferroni")
# significant between 6 degree and both other temps

# comparisons for time variable
pan2 %>%
  pairwise_t_test(
    o2_poc ~ day, paired = TRUE, 
    p.adjust.method = "bonferroni")
# significant between day 21 and day 27


## C:N ratio
# Check assumptions
pan2 <- pan %>% select(cn, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(cn, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "cn",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "cn", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = cn, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity not violated

res.aov_cn <- get_anova_table(res.aov)
# significant interaction + significant effects of day

# add main results to stat table
stat[7:9,3:6] <- res.aov_cn[,2:5]

# Considering a significant interaction
# main effect of day
one.way <- pan2 %>%
  group_by(temp) %>%
  anova_test(dv = cn, wid = plankto_ID, within = day) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# there is a significant main effect within 12 degree
# thus perform pairwise comparisons

# pairwise comparisons
pwc <- pan2 %>%
  group_by(temp) %>%
  pairwise_t_test(
    cn ~ day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)
# 6 degree differs between day 15 and 27 as well as between 24 and 27

## C:P ratio
# Check assumptions
pan2 <- pan %>% select(cp, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(cp, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "cp",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "cp", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = cp, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity not violated

res.aov_cp <- get_anova_table(res.aov)
# significant interaction + significant effects of day

# add main results to stat table
stat[10:12,3:6] <- res.aov_cp[,2:5]

# Considering a significant interaction
# main effect of day
one.way <- pan2 %>%
  group_by(temp) %>%
  anova_test(dv = cp, wid = plankto_ID, within = day) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# there is a significant main effect within all temperatures
# thus perform pairwise comparisons

# pairwise comparisons
pwc <- pan2 %>%
  group_by(temp) %>%
  pairwise_t_test(
    cp ~ day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)
# 12 degree differs between day 21 and 24 as well as between 24 and 27
# 18 degree differs between day 15 +18/21

write_xlsx(stat, "Data/StatisticsEFP.xlsx")

#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
require(NCmisc)
packages <- list.functions.in.file("~/AWI/RProjects/SOT22/EcosystemFunctions.R", alphabetic = TRUE) # set to your file path
summary(packages)
