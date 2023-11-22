### TopTrons mesozooplankton data ###
## Antonia Ahme, 22.11.2023 ##

### Load packages & housekeeping
require(dplyr)
require(tidyr)
require(vegan)
require(qualpalr)
require(ggplot2)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/TopTrons/Supplement")

### Create color palettes, themes & functions
class_pal <- qualpal(51, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

# Create the designs for plotting
bar.theme <- theme(panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title.x=element_blank(),
  axis.text.x = element_text(size = 20, face= "bold"),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  plot.title = element_blank(),
  strip.text = element_text(size = 30, face = "bold"),
  strip.background = element_blank(), strip.placement = "outside",
  text = element_text(size = 20, face = "bold"), legend.position = "right")

### Prepare the data
zoo <- readxl::read_excel("Data/MesoZoo.xlsx")
zoo$Temperature <- as.factor(zoo$Temperature)

### Bargraph on class level per treatment
# Rename classes for prettier plotting
zoo <- zoo %>%
  mutate(Species = recode(Species, "nauplii" = 'Nauplii')) %>%
  mutate(Species = recode(Species, "Asteroidea" = 'Asteroidea sp.')) %>%
  mutate(Species = recode(Species, "Bivalvia" = 'Bivaliva sp.')) %>%
  mutate(Species = recode(Species, "unidentified copepod" = 'Unidentified')) %>%
  mutate(Species = recode(Species, "unidentified zooplankton" = 'Unidentified')) %>%
  mutate(Species = recode(Species, "Pseudocalanus" = 'Pseudocalanus sp.')) %>%
  mutate(Species = recode(Species, "Acartia (Acartiura) clausi" = 'Acartia clausi')) %>%
  mutate(Species = recode(Species, "Paracalanus parvus parvus" = 'Paracalanus parvus'))

## Plotting
zoo_plot <- ggplot(zoo, aes(fill = Species, x = Replicate, y = Percent)) +
  facet_wrap(~ Temperature) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = class_pal$hex)

zoo_plot
ggsave("Output/MesoZoo.png", zoo_plot, height = 5, width = 9, dpi = 320)

### Bargraph of feeding categories per treatment
cat_plot <- ggplot(zoo, aes(fill = Food, x = Replicate, y = Percent)) +
  facet_wrap(~ Temperature) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = class_pal$hex)

cat_plot
ggsave("Output/MesoZooFeeding.png", cat_plot, height = 5, width = 9, dpi = 320)

#### Statistics on abundances ####
zoo_stat <- zoo %>%
  group_by(Temperature, Replicate) %>% 
  summarise(count=sum(Amount),
            .groups = 'drop')

kruskal.test(zoo_stat$count~zoo_stat$Temperature)  
