#### TOPTRONS REPLICATE ECOSYSTEM FUNCTIONS ####
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
require(ggrepel)

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
temp_pal <- c("P2"="skyblue2", "P5" = "steelblue2", "P6" = "mediumblue", "P10" = "royalblue1",
              "P3"="goldenrod2", "P8" = "goldenrod", "P9" = "yellow3", "P11" = "gold",
              "P1"="tomato2", "P4" = "salmon2", "P7" = "firebrick1", "P12" = "red")

# Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#### UPLOAD AND TIDY UP DATA ####
pan <- readxl::read_excel("Data/ParticulateNutrientsChlorophyll.xlsx")
pan$temp <- as.factor(pan$temp)
pan <- na.omit(pan)

oxy <- readxl::read_excel("Data/GrossOxygenProduction.xlsx")
oxy$temp <- as.factor(oxy$temp)
oxy <- subset(oxy, day >6)

din <- readxl::read_excel("Data/DissolvedNutrients.xlsx")
din$temp <- as.factor(din$temp)

dis <- readxl::read_excel("Data/DissolvedSilicate.xlsx")
dis$temp <- as.factor(dis$temp)

cc <- readxl::read_excel("Data/CarbonateChemistry.xlsx")
cc$temp <- as.factor(cc$temp)

pH <- readxl::read_excel("Data/MetaNoOut.xlsx")
pH$temp <- as.factor(pH$temp)


#### PLOT 1: CHLOROPHYLL A ####
# CV over time
chl <- ggplot(pan, aes(day, chl_µgL, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(pan, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("Chlorophyll a (µg/L)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3,30,3))

chl

ggsave("Output/ChlorophyllReps.png", chl, dpi = 300, width = 8, height = 4)

#### PLOT 2: ECOSYSTEM FUNCTIONS ####
# POC
p1 <- ggplot(pan, aes(day, µgL_c, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(pan, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("POC (µg/L)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3, 27, 3))

p1


# C:N
p2 <- ggplot(pan, aes(day, cn, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(pan, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("C:N ratio")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3, 27, 3))

p2

# C:P 
p3 <- ggplot(pan, aes(day, cp, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(pan, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("C:P ratio")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3, 27, 3))

p3

# GOP
p4 <- ggplot(oxy, aes(day, o2_poc, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(oxy, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("GOP")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(3, 30, 3))

p4

png("Output/EcosystemFunctionsReps.png", 4000, 4000, res = 400, type='cairo')
multiplot(p1, p2, p3, p4, cols=1)
dev.off()

#### PLOT 3: DISSOLVED NUTRIENTS ####
# Nitrate + nitrite
p1 <- ggplot(din, aes(day, no_µmolL, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(din, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("Nitrate + Nitrite (µM)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(0, 27, 3))

p1

# Phosphate
p2 <- ggplot(din, aes(day, po_µmolL, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(din, day == 9), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("Phosphate (µM)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(0, 27, 3))

p2

# Silicate
p3 <- ggplot(dis, aes(day, si_µM, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(dis, day == 28), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("Silicate (µM)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(0, 28, 2))

p3

# Export graphs as high resolution png file
png("Output/DissolvedNutrientsReps.png", 4000, 4000, res = 400, type='cairo')
multiplot(p1, p2, p3, cols=1)
dev.off()

#### PLOT 4: CARBONATE SYSTEM ####
# pH
p1 <- ggplot(pH, aes(day, pH_NBS, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(cc, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("pH (NBS)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(0, 27, 3))

p1

# DIC
p2 <- ggplot(cc, aes(day, DIC_mmolkgSW, col=plankto_ID)) + 
  geom_point(size = 3) + 
  geom_line(size = 0.75) +
  geom_label_repel(data = filter(cc, day == 27), 
                   aes(label = rep),
                   nudge_x = .75,
                   na.rm = TRUE) +
  plot.theme + theme(legend.position = "none",
                     axis.title.x = element_text(),
                     axis.title.y = element_text()) +
  labs(x="Incubation time (d)", y=bquote("DIC (mmol/kgSW)")) +
  scale_color_manual(values=temp_pal) +
  scale_x_continuous(breaks = seq(0, 27, 3))

p2

# Export graphs as high resolution png file
png("Output/CarbonateChemistryReps.png", 4000, 3000, res = 400, type='cairo')
multiplot(p1, p2, cols=1)
dev.off()
