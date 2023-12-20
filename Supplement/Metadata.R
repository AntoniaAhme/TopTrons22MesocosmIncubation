#### TOPTRONS METADATA ####
## Antonia Ahme, 11.01.2023 ##

rm(list=ls())
setwd("~/AWI/RProjects/TopTrons/Supplement")

#### LOAD PACKAGES ####
library(dplyr)
library(plyr)
library(readxl)
library(emmeans)
library(car)
library(nlme)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggpmisc)
library(ggplot2)

#### FUNCTIONS ####
# Functions for the data summary we to prepare data for graphs with mean and standard dev
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

#### UPLOAD AND TIDY UP DATA ####
meta <- readxl::read_excel("Data/Metadata.xlsx")
meta$temp <- as.factor(meta$temp)

# Create summaries of data
temp_deg <- data_summary(meta, varname="temp_DEG", 
                         groupnames=c("temp","day"))

sal_psu <- data_summary(meta, varname="sal_PSU", 
                        groupnames=c("temp","day"))

#### PLOT: TEMP, SAL####
# Plot real temperature over time for all temperatures
p1 <- ggplot(temp_deg, aes(x=day, y=temp_DEG, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=temp_DEG-sd, ymax=temp_DEG+sd), size=1.1, width=.1,
                position=position_dodge(0.05)) +
  plot.theme+
  labs(color = "Temperature (°C)") +
  labs(x="Incubation time (d)", y=bquote("Temperature (°C)")) + 
  scale_color_manual(values=c("6" = "skyblue2", "12" = "goldenrod1" , "18" = "tomato2")) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 20, 2))

p1

ggsave("Output/Temperature.png", p1, dpi = 300, width = 8, height = 4)

# Plot salinity over time for all temperatures
p2 <- ggplot(sal_psu, aes(x=day, y=sal_PSU, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=sal_PSU-sd, ymax=sal_PSU+sd), size=1.1, width=.1,
                position=position_dodge(0.05)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position="none") +
  labs(color = "Temperature [°C]") +
  labs(x="Incubation time (d)", y=bquote("Salinity (PSU)")) + 
  scale_color_manual(values=c("6" = "skyblue2", "12" = "goldenrod1" , "18" = "tomato2")) +
  scale_y_continuous(limits = c(25,35))

p2

ggsave("Output/Salinity.png", p2, dpi = 300, width = 8, height = 4)
