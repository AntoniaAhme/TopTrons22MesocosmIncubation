### TOPTRONS CHLOROPHYLL NUTRIENTS ###
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

#### UPLOAD AND TIDY UP DATA ####
chl <- readxl::read_excel("Data/Chlorophyll.xlsx")
chl$temp <- as.factor(chl$temp)
chl <- subset(chl, day > 12)
chl_full <- chl

din <- readxl::read_excel("Data/DissolvedNutrients.xlsx")
din$temp <- as.factor(din$temp)
din <- subset(din, day > 12)
din_full <- din

dis <- readxl::read_excel("Data/DissolvedSilicate.xlsx")
dis$temp <- as.factor(dis$temp)
dis <- subset(dis, day > 12)
dis_full <- dis

#### PLOTTING ####
# Chlorophyll
chl <- data_summary(chl, varname="chl_µgL", 
                    groupnames=c("temp","day"))

p1 <- ggplot(chl, aes(x=day, y=chl_µgL, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=chl_µgL-sd, ymax=chl_µgL+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('Chlorophyll' ~ '(' ~ 'µg L' ^ -1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 16, 2))+
  scale_color_manual(values=temp_pal)

p1

# Nitrate
din <- data_summary(din_full, varname="no_µmolL", 
                    groupnames=c("temp","day"))

p2 <- ggplot(din, aes(x=day, y=no_µmolL, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=no_µmolL-sd, ymax=no_µmolL+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('Nitrate' ~ '(' ~ 'µmol L' ^ -1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 20, 5))+
  scale_color_manual(values=temp_pal)

p2


# Phosphate
dip <- data_summary(din_full, varname="po_µmolL", 
                    groupnames=c("temp","day"))

p3 <- ggplot(dip, aes(x=day, y=po_µmolL, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=po_µmolL-sd, ymax=po_µmolL+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('Phospate' ~ '(' ~ 'µmol L' ^ -1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 0.3, 0.05))+
  scale_color_manual(values=temp_pal)

p3

# Phosphate
dis <- data_summary(dis, varname="si_µM", 
                    groupnames=c("temp","day"))

p4 <- ggplot(dis, aes(x=day, y=si_µM, color=temp, group=temp)) + 
  geom_point(position=position_dodge(0.5), size = 4)+
  geom_line(position=position_dodge(0.5), size = 1)+
  geom_errorbar(aes(ymin=si_µM-sd, ymax=si_µM+sd), size=1.1, width=.8,
                position=position_dodge(0.5)) +
  plot.theme +
  labs(x=bquote('Incubation time (d)'), y=bquote('Silicate' ~ '(' ~ 'µmol L' ^ -1 ~ ')')) +
  scale_x_continuous(breaks = seq(0, 27, 3))+
  scale_y_continuous(breaks = seq(0, 20, 2))+
  scale_color_manual(values=temp_pal)

p4

# Export graphs as high resolution png file
png("Output/ChlorophyllNutrients.png", 5000, 4000, res = 400, type='cairo')
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

#### STATISTICS ####
### TWO-WAY REPEATED MEASURES ANOVA
#https://stats.stackexchange.com/questions/181563/analyzing-repeated-measures-experiment-with-multiple-treatment-groups-and-multip
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
#https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

## Create a table for statistical outout of two-way RM anova
Parameter <- rep(c("Chlorophyll", "Nitrate", "Phosphate", "Silicate"), each=3)
Effect <- rep(c("Temperature","Time","Temperature:Time"), times=4)
stat <- data.frame(Parameter, Effect)
stat['DFn'] <- NA
stat['DFd'] <- NA
stat['F'] <- NA
stat['P'] <- NA

## Chlorophyll
# Check assumptions
pan2 <- chl_full %>% select(chl_µgL, temp, day, plankto_ID)
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

res.aov_chl <- get_anova_table(res.aov)
# no significant interaction + significant effects of day and temp

# add main results to stat table
stat[1:3,3:6] <- res.aov_chl[,2:5]

## Nitrate
# Check assumptions
pan2 <- din_full %>% select(no_µmolL, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(no_µmolL, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "no_µmolL",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "no_µmolL", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = no_µmolL, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity violated

res.aov_din <- get_anova_table(res.aov)
# no significant interaction + significant effects of day

# add main results to stat table
stat[4:6,3:6] <- res.aov_din[,2:5]

## Phosphate
# Check assumptions
pan2 <- din_full %>% select(po_µmolL, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(po_µmolL, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "po_µmolL",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "po_µmolL", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but looks all right

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = po_µmolL, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# assumption of sphericity not met - GG correction applied

res.aov_dip <- get_anova_table(res.aov)
# no significant interaction + significant effects of day

# add main results to stat table
stat[7:9,3:6] <- res.aov_dip[,2:5]

# Considering no significant interaction
# comparisons for time variable
pan2 %>%
  pairwise_t_test(
    po_µmolL ~ day, paired = TRUE, 
    p.adjust.method = "bonferroni")

## Silicate
# Check assumptions
pan2 <- dis_full %>% select(si_µM, temp, day, plankto_ID)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(si_µM, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "si_µM",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "si_µM", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = si_µM, wid = plankto_ID, 
  between = temp, within = day)
res.aov
# sphericity not violated

res.aov_si <- get_anova_table(res.aov)
# significant interaction + significant effects of day and temp

# add main results to stat table
stat[10:12,3:6] <- res.aov_si[,2:5]

# Considering a significant interaction
# main effect of temp
one.way <- pan2 %>%
  group_by(day) %>%
  anova_test(dv = si_µM, wid = plankto_ID, between = temp) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# there is a significant main effect for temp on days 20-28
# thus perform pairwise comparisons

# pairwise comparisons
pwc <- pan2 %>%
  group_by(day) %>%
  pairwise_t_test(
    si_µM ~ temp, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)

# main effect of day
one.way <- pan2 %>%
  group_by(temp) %>%
  anova_test(dv = si_µM, wid = plankto_ID, within = day) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# 
# thus perform pairwise comparisons

# pairwise comparisons
pwc <- pan2 %>%
  group_by(temp) %>%
  pairwise_t_test(
    si_µM ~ day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)

## Save main results from rmANOVA
write_xlsx(stat, "Data/StatisticsChlNut.xlsx")
