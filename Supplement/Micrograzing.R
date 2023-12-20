#### TopTrons Micrograzing ####
## Antonia Ahme, 27.03.2023 ##

rm(list=ls())
setwd("~/AWI/RProjects/TopTrons/Supplement")

#### LOAD PACKAGES ####
require(dplyr)
require(plyr)
require(readxl)
require(emmeans)
require(car)
require(nlme)
require(tidyverse)
require(ggpubr)
require(rstatix)
require(ggpmisc)
require(ggplot2)
require(writexl)

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

# CV summary function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#### UPLOAD AND TIDY UP DATA ####
gra <- readxl::read_excel("Data/GrazingRates.xlsx")
gra$Temperature <- as.factor(gra$Temperature)

#### SUMMARY DATASETS FOR PLOTTING ####
gra_sum <- data_summary(gra, varname="m_d", 
                        groupnames=c("Temperature","Day"))


#### PLOT 1: GRAZING RATES ####
gra_plot <- ggplot(gra_sum, aes(x=Day, y=m_d, color = Temperature, group=Temperature)) + 
  geom_point(position=position_dodge(2), size = 3)+
  geom_errorbar(aes(ymin=m_d-sd, ymax=m_d+sd), size=1.1, width=1.25,
                position=position_dodge(2)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position="none") +
  labs(x="Incubation time (d)", y=bquote('m' ~ '(' ~ d ^ -1 ~ ')')) + 
  scale_x_continuous(breaks = seq(15, 27, 12))+
  scale_color_manual(values=c("6" = "skyblue2", "12" = "goldenrod1" , "18" = "tomato2"))

gra_plot

ggsave("Output/Micrograzing.png", gra_plot, dpi = 300, width = 6, height = 5)

#### STATISTICS ####
### TWO-WAY REPEATED MEASURES ANOVA
#https://stats.stackexchange.com/questions/181563/analyzing-repeated-measures-experiment-with-multiple-treatment-groups-and-multip
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
#https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

## Shannon
# Check assumptions
pan2 <- gra %>% select(m_d, Temperature, Day, Replicate)
pan2$Day <- as.factor(pan2$Day)
pan2$Replicate <- as.factor(pan2$Replicate)

# Summary
pan2 %>%
  group_by(Temperature, Day) %>%
  get_summary_stats(m_d, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "Day", y = "m_d",
  color = "Temperature", palette = temp_pal)
bxp

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = m_d, wid = Replicate, 
  between = Temperature, within = Day)
res.aov
# sphericity not violated and only interaction of time and temp

# pairwise comparisons of time
pwc <- pan2 %>%
  group_by(Temperature) %>%
  pairwise_t_test(
    m_d ~ Day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)

# pairwise comparisons of temperatures
pwc <- pan2 %>%
  group_by(Day) %>%
  pairwise_t_test(
    m_d ~ Temperature, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)
