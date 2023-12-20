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
                    axis.title.x = element_text(size = 15),
                    axis.text.x = element_text(size = 15, face= "bold"),
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_text(size = 15),
                    axis.text.y = element_text(size = 15, face= "bold"),
                    axis.ticks.y = element_blank(),
                    plot.title = element_text(size = 15, face= "bold"),
                    strip.text = element_text(size = 15, face = "bold"),
                    strip.background = element_blank(), strip.placement = "outside",
                    text = element_text(size = 15, face = "bold"), legend.position = "none")

# Function for the data summary we to prepare data for graphs with mean and sd
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

# Color palette
temp_pal <- c("6"="skyblue2", "12"="goldenrod2", "18"="tomato2")
#temp_pal <- c("6"="darkorange", "12"="blue", "18"="turquoise3") # inverted for presentation
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
poc <- data_summary(pan, varname="µgL_c", 
                    groupnames=c("temp","day"))

p1 <- ggplot(poc, aes(x=day, y=µgL_c, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=µgL_c-sd, ymax=µgL_c+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('POC' ~ '(' ~ 'µg L' ^ -1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 4000, 500))+
  scale_color_manual(values=temp_pal)

p1

# GOP
gop <- data_summary(oxy, varname="o2_poc", 
                    groupnames=c("temp","day"))

p2 <- ggplot(gop, aes(x=day, y=o2_poc, color=temp, group=temp,day)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=o2_poc-sd, ymax=o2_poc+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('GOP' ~ '(' ~ 'µmol O' [2] ~ 'mg POC' ^ -1 ~ 'd' ^-1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 30, 5))+
  scale_color_manual(values=temp_pal)

p2

# C:N ratio
c_n <- data_summary(pan, varname="cn", 
                    groupnames=c("temp","day"))

p3 <- ggplot(c_n, aes(x=day, y=cn, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=cn-sd, ymax=cn+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote('C:N ratio')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 45, 5))+
  scale_color_manual(values=temp_pal)

p3

# C:P ratio
c_p <- data_summary(pan, varname="cp", 
                    groupnames=c("temp","day"))

p4 <- ggplot(c_p, aes(x=day, y=cp, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote('C:P ratio')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 1000, 200))+
  scale_color_manual(values=temp_pal)

p4

# Export graphs as high resolution png file
png("Output/EcosystemFunctions.png", 5000, 4000, res = 400, type='cairo')
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

#### STATISTICS ####
### TWO-WAY REPEATED MEASURES ANOVA
#https://stats.stackexchange.com/questions/181563/analyzing-repeated-measures-experiment-with-multiple-treatment-groups-and-multip
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
#https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

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
#ggqqplot(pan2, "µgL_c", ggtheme = theme_bw()) +
#  facet_grid(day ~ temp, labeller = "label_both")
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
#ggqqplot(pan2, "o2_poc", ggtheme = theme_bw()) +
#  facet_grid(day ~ temp, labeller = "label_both")
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
#ggqqplot(pan2, "cn", ggtheme = theme_bw()) +
# facet_grid(day ~ temp, labeller = "label_both")
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
#ggqqplot(pan2, "cp", ggtheme = theme_bw()) +
#  facet_grid(day ~ temp, labeller = "label_both")
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

write_xlsx(stat, "Data/StatisticsEcoFunctions.xlsx")

#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
require(NCmisc)
packages <- list.functions.in.file("~/AWI/RProjects/SOT22/EcosystemFunctions.R", alphabetic = TRUE) # set to your file path
summary(packages)
