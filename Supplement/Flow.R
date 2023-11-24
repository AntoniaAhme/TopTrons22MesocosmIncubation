#### TopTrons Flow Cytometry Data - Phycocyanin-containing bacteria ####
## Antonia Ahme, 21.11.2023 ##

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
require(qualpalr)
require(ggrepel)

#### LAYOUTS ####
# Create plot theme
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
temp_pal <- c("P2"="skyblue2", "P5" = "steelblue2", "P6" = "mediumblue", "P10" = "royalblue1",
              "P3"="goldenrod2", "P8" = "goldenrod", "P9" = "yellow3", "P11" = "gold",
              "P1"="tomato2", "P4" = "salmon2", "P7" = "firebrick1", "P12" = "red")

#### UPLOAD AND TIDY UP DATA ####
count <- readxl::read_excel("Data/FlowCounts.xlsx")
count$temp <- as.factor(count$temp)
count <- na.omit(count)
count <- count[, c("sample_ID", "plankto_ID", "temp", "rep", "day", "S_HFL-4")]
count <- plyr::rename(count,  c("S_HFL-4" = 'Abundance'))


#### PLOT BARPLOT OF PSEUDOMONAS PER MESOCOSM OVER TIME ####
phyco <- ggplot(count, aes(day, Abundance, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(count, day == 21), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("Abundance (mL^-1)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3, 27, 3))

phyco
ggsave("Output/Phycocyanin.png", phyco, height = 5, width = 9, dpi = 320)
