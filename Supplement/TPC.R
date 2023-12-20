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

# Plot theme
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
#ggsave("Output/TPC.png", tpc, width = 15, height =15 , dpi = 500, unit = "cm", device = "png")
#_______________________________________________________________________________

# Non-parametric Bootstrapping

library(rTPC)
library(nls.multstart)
library(broom)
library(gridExtra)
library(boot)
library(car)
library(patchwork)
library(minpack.lm)

#Bootstrap Analysis for CIs
fit_tpc_boot <- minpack.lm::nlsLM(r~thomas_2017(temp = Temp, a,b,c,d,e),
                                 data = df_fitting,
                                 start = start_vals,
                                 lower = low_lims,
                                 upper = upper_lims,
                                 weights = rep(1, times = nrow(df_fitting)))

# bootstrap using case resampling
boot1 <- Boot(fit_tpc_boot, method = 'case')

# look at the data
head(boot1$t)

boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(df_fitting$Temp), max(df_fitting$Temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = thomas_2017(temp, a,b,c,d,e))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()
boot1_conf_preds <-boot1_conf_preds
#

# calculate parameter CIs
extra_params <- calc_params(fit_tpc_boot) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_tpc_boot, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_H5_boot)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params <- left_join(ci_extra_params, extra_params)
#ci_extra_params <- mutate(ci_extra_params_H5, strain = "H5")


tpc_points <- read_xlsx("TPC_points.xlsx")
tpc_points$Temp <- as.numeric(tpc_points$Temp)


p1 <- ggplot() +
  geom_hline(yintercept = 0,linetype ="dashed")+
  geom_point(aes(Temp, r), tpc_points, size = 3, alpha = 1) +
  geom_line(aes(Temp, .fitted), preds, size = 1.5, show.legend = FALSE, color="black")+
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, alpha = 0.3)+
  scale_x_continuous(breaks = seq(0,30, by =5))+
  plot.theme+
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate')

p1
ggsave("Output/TPC.png", p1, width = 15, height =15 , dpi = 500, unit = "cm", device = "png")
