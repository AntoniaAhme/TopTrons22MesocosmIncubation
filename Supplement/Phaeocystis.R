### TOPTRONS PHAEOCYSTIS CORERLATIONS ###
## Antonia Ahme, 06.01.2024 ##

#### HOUSEKEEPING ####
require(dplyr)
require(vegan)
require(ggplot2)
require(ggpubr)
require(ggforce)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/TopTrons/Supplement") # set to your working directory

set.seed(22)

#### CREATE LAYOUTS & FUNCTIONS ####
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

# Color palettes
temp_pal <- c("6" = "skyblue2", "12" = "goldenrod2", "18" = "tomato2")

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

#### UPLOAD DATA ####
## Load Phaeo dataframe
phaeo <- readxl::read_excel("Data/Phaeocystis.xlsx")

## Load variables
pan <- readxl::read_excel("Data/ParticulateNutrients.xlsx")
pan$temp <- as.factor(pan$temp)
pan <- subset(pan, day > 12)
pan$Sample <- paste(pan$plankto_ID, pan$tp, sep = "-")
pan <- pan %>% select(Sample, µgL_c, cp)

din <- readxl::read_excel("Data/DissolvedNutrients.xlsx")
din$temp <- as.factor(din$temp)
din <- subset(din, day > 12)
din$Sample <- paste(din$plankto_ID, din$tp, sep = "-")
din <- din %>% select(Sample, no_µmolL)

cc <- readxl::read_excel("Data/Metadata.xlsx")
cc$temp <- as.factor(cc$temp)
cc <- subset(cc, day > 12)
cc$Sample <- paste(cc$plankto_ID, cc$tp, sep = "-")
cc <- cc %>% select(Sample, pH_NBS)

## Add variables to phaeo dataframe
phaeo_all <- subset(phaeo, day > 12)
phaeo_all <- merge(phaeo_all, pan, by = "Sample")
phaeo_all <- merge(phaeo_all, din, by = "Sample")
phaeo_all <- merge(phaeo_all, cc, by = "Sample", all = TRUE)

phaeo_all$day <- as.factor(phaeo_all$day)

#### PLOTTING ####
## Biomass
p1 <- ggplot(phaeo_all, aes(x=Abundance, y=µgL_c)) + 
  geom_point(position=position_dodge(0.5), size = 4, aes(color=temp, shape=day))+
  plot.theme +
  geom_smooth(method=lm, color = "black") +
  stat_regline_equation(label.x = 3, label.y = 3000)+
  stat_cor(label.x = 3, label.y = 3150) +
  labs(x=bquote('Normalised read abundance of P. globosa'), y=bquote('POC' ~ '(' ~ 'µg L' ^ -1 ~ ')')) +
  scale_color_manual(values=temp_pal)

p1

## C:P
p2 <- ggplot(phaeo_all, aes(x=Abundance, y=cp)) + 
  geom_point(position=position_dodge(0.5), size = 4, aes(color=temp, shape=day))+
  plot.theme +
  geom_smooth(method=lm, color = "black") +
  stat_regline_equation(label.x = 3, label.y = 950)+
  stat_cor(label.x = 3, label.y = 1000) +
  labs(x=bquote('Normalised read abundance of P. globosa'), y=bquote('C:P')) +
  scale_color_manual(values=temp_pal)

p2

## Nitrate
p3 <- ggplot(phaeo_all, aes(x=Abundance, y=no_µmolL)) + 
  geom_point(position=position_dodge(0.5), size = 4, aes(color=temp, shape=day))+
  plot.theme +
  geom_smooth(method=lm, color = "black") +
  stat_regline_equation(label.x = 3, label.y = 1)+
  stat_cor(label.x = 3, label.y = 2.5) +
  labs(x=bquote('Normalised read abundance of P. globosa'), y=bquote('Nitrate' ~ '(' ~ 'µmol L' ^ -1 ~ ')')) +
  scale_color_manual(values=temp_pal)

p3

## pH
p4 <- ggplot(phaeo_all, aes(x=Abundance, y=pH_NBS)) + 
  geom_point(position=position_dodge(0.5), size = 4, aes(color=temp, shape=day))+
  plot.theme +
  geom_smooth(method=lm, color = "black") +
  stat_regline_equation(label.x = 3, label.y = 8)+
  stat_cor(label.x = 3, label.y = 8.03) +
  labs(x=bquote('Normalised read abundance of P. globosa'), y=bquote('pH (NBS)')) +
  scale_color_manual(values=temp_pal)

p4

#Export graphs as high resolution png file
png("Output/Phaeocystis.png", 4500, 3500, res = 400, type='cairo')
multiplot(p1, p2, p3, p4, cols=2)
dev.off()
