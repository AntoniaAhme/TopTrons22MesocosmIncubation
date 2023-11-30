### TOPTRONS TPC ###
## Jakob Giesler, 30.11.2023 ##

#_______________________________________________________________________________

# Housekeeping
rm(list = ls())

setwd("~/AWI/RProjects/TopTrons/Supplement") # set to your working directory

library(tidyverse)
library(readxl)
library(growthrates)
library(ggthemes)
library(writexl)

#_______________________________________________________________________________

#### RUN 1 ####

# Data import & edit
tpc_data <- read_excel("Data/TPCFluorescence.xlsx")
tpc_data <- select(tpc_data, ID_Unique, Temp, Timepassed, Chl)
tpc_data <- mutate(tpc_data, r = NA, p_value = NA) #Write empty columns to fill
tpc_data$Chl <- tpc_data$Chl+1
tpc_data <- group_by(tpc_data, Temp, Timepassed)
tpc_data <- mutate(tpc_data, mean_chl = mean(Chl))
tpc_data$Temp <- as.factor(tpc_data$Temp)

# First look at Chl- fluorescence over time
mean_tpc <- tpc_data %>% ggplot(aes(Timepassed, mean_chl, color = Temp))+
  geom_point()+
  geom_line()+
  theme_classic()+
  scale_color_viridis_d()+
  ylab("mean Chl (FU)")+
  xlab("Time (d)")
mean_tpc

str(tpc_data)

tpc_data <- select(tpc_data, ID_Unique, Temp, Chl, Timepassed)
tpc_data <- filter(tpc_data, Timepassed == "1" | Timepassed == "5")
tpc_data <- spread(tpc_data, key = Timepassed, value = Chl)
tpc_data <- rename(tpc_data, "Start" = "1")
tpc_data <- rename(tpc_data, "End" = "5")
tpc_data <- mutate(tpc_data, r = (log(End)-log(Start))/4)

tpc_data %>% ggplot(aes(Temp, r, group = as.factor(Temp)))+
  geom_point()+
  theme_classic()

# safe points for later to add the real data points to the modeled TPCs
write_xlsx(tpc_data, "Data/TPC_points.xlsx")

#_______________________________________________________________________________

# TPC fitting

library(rTPC)
library(nls.multstart)
library(broom)
library(gridExtra)

# prepare data for TPC fitting
df_fitting <- tpc_data
names(df_fitting)
df_fitting <- dplyr::ungroup(df_fitting)
df_fitting <- dplyr::select(df_fitting,Temp, ID_Unique, r)
df_fitting <- dplyr::group_by(df_fitting, Temp)
df_fitting <- mutate(df_fitting, curve_id = 1)

# check points
df_fitting %>% ggplot(aes(Temp, r))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'µ (day^-1)',
       title = 'µ across temperatures')

df_fitting$Temp <- as.numeric(df_fitting$Temp)
df_fitting$Temp <- df_fitting$Temp*3

# choose model
mod = 'thomas_2017'

# get start values
start_vals <- get_start_vals(df_fitting$Temp, df_fitting$r, model_name = 'thomas_2017')

# get limits
low_lims <- get_lower_lims(df_fitting$Temp, df_fitting$r, model_name = 'thomas_2017')
upper_lims <- get_upper_lims(df_fitting$Temp, df_fitting$r, model_name = 'thomas_2017')

start_vals
low_lims
upper_lims

# fit model
fit <- nls_multstart(r~thomas_2017(temp = Temp, a,b,c,d,e),
                        data = df_fitting,
                        iter = c(3,3,3,3,3),
                        start_lower = start_vals - 10,
                        start_upper = start_vals + 10,
                        lower = get_lower_lims(df_fitting$Temp, df_fitting$r, model_name = 'thomas_2017'),
                        upper = get_upper_lims(df_fitting$Temp, df_fitting$r, model_name = 'thomas_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)
summary(fit)

# get predictions
preds <- data.frame(Temp = seq(min(df_fitting$Temp), max(df_fitting$Temp), length.out = 100))
preds <- broom::augment(fit, newdata = preds)

# plot TPC
tpc <- ggplot(preds) +
  geom_point(aes(Temp, r), df_fitting) +
  geom_line(aes(Temp, .fitted), col = 'blue') +
  theme_few()+
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 3)) +
  ylab(bquote('µ' ~ '(' ~ d ^ -1 ~ ')')) +
  xlab('Temperature (°C)')+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
tpc
ggsave("Output/TPC.png", tpc, width = 15, height =15 , dpi = 500, unit = "cm", device = "png")
#_______________________________________________________________________________



