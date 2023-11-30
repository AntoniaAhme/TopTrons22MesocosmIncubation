### TOPTRONS COMMUNITY COMPOSITION ###
## Antonia Ahme, 14.11.2023 ##

#### HOUSEKEEPING ####
require(dplyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(zCompositions)
require(propr)
require(ggpubr)
require(rstatix)
require(microbial)
require(qualpalr)
require(microViz)
require(ggforce)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/TopTrons") # set to your working directory

set.seed(22)

#### CREATE LAYOUTS & FUNCTIONS ####
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
  strip.text = element_text(size = 15, face = "bold"),
  strip.background = element_blank(), strip.placement = "outside",
  text = element_text(size = 20, face = "bold"), legend.position = "right")

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
temp_pal <- c("skyblue2", "goldenrod2", "tomato2")
temp_pal2 <- c("lightskyblue1", "gold", "salmon1")
uni_pal <- c("gold", "gold", "gold", "gold", "gold", "salmon1", "salmon1",
             "salmon1",  "salmon1", "salmon1", "lightskyblue1", "lightskyblue1",
              "lightskyblue1", "lightskyblue1", "lightskyblue1")

# Function for italicizing genus legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

# Standard error and mean function
se <- function(x, ...) sqrt(var(x, ...)/length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#### UPLOAD DATA ####
# Import asv and tax tab
asv <- read.delim('Data/Counts.txt')
asv <- data.frame(asv, row.names = 1)
colnames(asv) <- sub("SOT22.", "", colnames(asv))
colnames(asv) <- sub(".euk", "", colnames(asv))
colnames(asv) <- gsub('\\.', '-', colnames(asv))
tax <- read.delim('Data/Taxonomy.txt')

# Import sample tab
sam <- read.delim('Data/Samples.txt')
sam$sample_ID <- sub("SOT22-", "", sam$sample_ID)
sam$sample_ID <- sub("-euk", "", sam$sample_ID)
sam$sample_ID <- as.factor(sam$sample_ID)
sam$day <- as.numeric(sam$day)
sam$temp <- as.factor(sam$temp)

#### CLEAN, NORMALISE & TRANSFORM DATA ####
asv_all <- asv
sam_all <- sam

# Remove ASVs that don't occur in the dataset
asv <- asv[rowSums(asv)>0,]

# Investigate sequencing depth and remove samples where it is too low
depths <- colSums(asv) # prepare df containing the sequencing depth
plot(depths) # sequencing depth does not have any outliers or very low numbers
min(depths)
max(depths)
quantile(depths, probs = seq(0,1,.5)) # Show all sequencing depth quantiles
# Remove all samples outside of the 90% quantile range
asv <- asv[,colSums(asv)<272993]
asv <- asv[,colSums(asv)>32233]
rm(depths)

# Match sam tab to new asv tab as some samples might have been removed in the last step
sam <- sam[sam$sample_ID %in% colnames(asv),]

# Check and adjust whether rownames of meta-info file match colnames of ASV counts
all(sam$sample_ID == colnames(asv))

# Remove asvs with a count of less than 10 (singletons) in replicate sample means (create rep means first)
ASV <- asv
colnames(ASV) <- sam$uni_ID[match(colnames(ASV),sam$sample_ID)]
Names <- unique(names(ASV))
ASV <- sapply(Names, function(x)  rowMeans(ASV[names(ASV) %in% x]))
ASV <- as.data.frame(ASV)
ASV <- ASV %>% filter_at(vars(1:28), any_vars(.>=10)) # number of pooled samples and cutoff-level of 10
ASV.rn <- rownames(ASV)
asv <- asv[rownames(asv) %in% ASV.rn,]
rm(ASV,Names,ASV.rn)

tax <- tax[tax$ASV %in% rownames(asv),] # create matching taxonomy table after removal

# Check taxonomy
unique(tax$Supergroup)
unique(tax$Division)
unique(tax$Class)

# Replace taxonomic levels that were not annotated with "other"
tax[is.na(tax)] <- "Other"

# Remove unwanted groups
tax <- filter(tax, tax$Division!="Metazoa")
tax <- filter(tax, tax$Division!="Fungi")
tax <- filter(tax, tax$Division!="Cryptophyta:nucl")
tax <- filter(tax, tax$Division!="Chlorophyta:plas")
tax <- filter(tax, tax$Division!="Ochrophyta:plas")
tax <- filter(tax, tax$Supergroup!="Archaeplastida:plas")
tax <- filter(tax, tax$Supergroup!="Hacrobia:nucl")

tax_all <- tax

# Create tax file for phytoplankton by exluding hetero- and mixotrophs
# based on https://doi.org/10.1111/jeu.12691
tax <- filter(tax, tax$Supergroup!="Cercozoa")
tax <- filter(tax, tax$Supergroup!="Apusozoa")
tax <- filter(tax, tax$Supergroup!="Amoebozoa")
tax <- filter(tax, tax$Supergroup!="Rhizaria")
tax <- filter(tax, tax$Division!="Choanoflagellida")
tax <- filter(tax, tax$Division!="Pseudofungi")
tax <- filter(tax, tax$Division!="Opalozoa")
tax <- filter(tax, tax$Division!="Picozoa")
tax <- filter(tax, tax$Division!="Telonemia")
tax <- filter(tax, tax$Division!="Katablepharidophyta")
tax <- filter(tax, tax$Division!="Sagenista")
tax <- filter(tax, tax$Division!="Ciliophora")
tax <- filter(tax, tax$Division!="Apicomplexa")
tax <- filter(tax, tax$Division!="Mesomycetozoa")
tax <- filter(tax, tax$Division!="Centroheliozoa")
tax <- filter(tax, tax$Division!="Stramenopiles_X")
tax <- filter(tax, tax$Division!="Streptophyta")
tax <- filter(tax, tax$Division!="Dinoflagellata") 
tax <- filter(tax, tax$Class!="Syndiniales")
tax <- filter(tax, tax$Class!="Noctilucophyceae")
tax <- filter(tax, tax$Class!="Chrysophyceae")

unique(tax$Division)
unique(tax$Class)

# Subset asv tab based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax$ASV,]

# Create a tax file for raw read phyloseq object
tax_raw <- tax

## Create a phyloseq object of the raw data for rarefaction curves
rownames(tax_raw) <- tax_raw$ASV
tax_raw <- tax_raw[,2:9] # delete ASV column
tax_raw <- as.matrix(tax_raw)
rownames(sam) <- sam$sample_ID
OTU = otu_table(asv, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_raw)
sam$uni_ID <- as.factor(sam$uni_ID)
samples = sample_data(sam)

ps_raw <- phyloseq(OTU, TAX, samples)

# Check rarefaction curve for justifying sampling depths
rarecurve(t(otu_table(ps_raw)), step=50, cex=0.5) # looks fine
# save manually in RSTudio

## Scaling with ranked subsampling (srs)
# All samples will be scaled to sample with lowest sequencing depth
depth.min <- min(colSums(asv))
asv.srs <- SRS(asv, depth.min) # running the SRS function
rownames(asv.srs) <- rownames(asv)
rm(depth.min)

# Check and adjust that rownames sam tab and colnames of asv tab match
all(rownames(sam) == colnames(asv))

# Prepare tax file for phyloseq object
rownames(tax) <- tax$ASV
tax <- tax[,2:9] # delete ASV column
tax <- as.matrix(tax)

# Name elements for phyloseq object with scaled data
OTU = otu_table(asv.srs, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)

# Create phyloseq object
ps <- phyloseq(OTU, TAX, SAM)

## Create CLR transformed data for beta diversity analyses
# Remove zeroes
czm <- cmultRepl(t(asv),  label=0, method="CZM") 

# Clr transformation
rho <- propr(czm, metric = 'rho', ivar = 'clr', symmetrize = TRUE, p=0)

# Clr data
clr <- rho@logratio

# Transpose the normalised and filtered asv tab & create dataframe for r
asv.clr <- t(clr)
asv.clr <- as.data.frame(asv.clr)
OTU = otu_table(asv.clr, taxa_are_rows = TRUE)

# Create a phyloseq object
ps_clr <- phyloseq(OTU, TAX, SAM)

#### CALCULATE & PLOT DIVERSITY ####
ps.rich <- microbial::richness(ps,  method = c("Observed", "Evenness", "Shannon"))

# Add diversity measures to sample tab
sam$Richness <- ps.rich$Observed
sam$Evenness <- ps.rich$Evenness
sam$Shannon <- ps.rich$Shannon
sam$cosm <- as.factor(sam$cosm)
div <- subset(sam, day > 12)

## Plot diversity
# Create summary of data for shannon index
div_shan <- data_summary(div, varname="Shannon", 
                         groupnames=c("day", "temp"))

shan_time <- ggplot(div_shan, aes(x=day, y=Shannon, color=temp)) + 
  geom_point(position=position_dodge(0.05), size = 4)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se), size=1.1, width=.8,
                position=position_dodge(0.05)) +
  plot.theme +
  labs(x="Incubation day", y=bquote("Shannon index")) + 
  scale_x_continuous(breaks = seq(15, 27, 3))+
  scale_color_manual(values=temp_pal)

shan_time
ggsave("Output/ShannonOverTime.png", shan_time, height = 5, width = 8, dpi = 320)

# Richness
# Create summary of data for species richness
div_rich <- data_summary(div, varname="Richness", 
                         groupnames=c("day", "temp"))

rich_time <- ggplot(div_rich, aes(x=day, y=Richness, color=temp)) + 
  geom_point(position=position_dodge(0.05), size = 4)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), size=1.1, width=.8,
                position=position_dodge(0.05)) +
  plot.theme +
  scale_x_continuous(breaks = seq(15, 27, 3))+
  labs(x="Day", y=bquote("Species Richness")) + 
  scale_color_manual(values=temp_pal)

rich_time
ggsave("Output/RichnessOverTime.png", rich_time, height = 5, width = 8, dpi = 320)

# Evenness
# Create summary of data for species evenness
div_even <- data_summary(div, varname="Evenness", 
                         groupnames=c("day", "temp"))

even_time <- ggplot(div_even, aes(x=day, y=Evenness, color=temp)) + 
  geom_point(position=position_dodge(0.05), size = 4)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Evenness-se, ymax=Evenness+se), size=1.1, width=.8,
                position=position_dodge(0.05)) +
  plot.theme +
  scale_x_continuous(breaks = seq(15, 27, 3))+
  labs(x="Day", y=bquote("Species Evenness")) + 
  scale_color_manual(values=temp_pal)

even_time
ggsave("Output/EvennessOverTime.png", even_time, height = 5, width = 8, dpi = 320)


### Statistics of diversity metrices
### TWO-WAY REPEATED MEASURES ANOVA
#https://stats.stackexchange.com/questions/181563/analyzing-repeated-measures-experiment-with-multiple-treatment-groups-and-multip
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
#https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

# Exclude groups with too few datapoints
div <- subset(div, day != "24")

## Shannon
# Check assumptions
pan2 <- div %>% select(Shannon, temp, day, cosm)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(Shannon, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "Shannon",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "Shannon", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = Shannon, wid = cosm, 
  between = temp, within = day)
res.aov
# sphericity not violated

# Considering a significant interaction
# main effect of day
one.way <- pan2 %>%
  group_by(temp) %>%
  anova_test(dv = Shannon, wid = cosm, within = day) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# there is a significant main effect within 6 degree
# thus perform pairwise comparisons

# pairwise comparisons
pwc <- pan2 %>%
  group_by(temp) %>%
  pairwise_t_test(
    Shannon ~ day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
print(pwc, n=30)
# 6 degree differs between day 15 and 27 as well as between 24 and 27


## Richness
# Check assumptions
pan2 <- div %>% select(Richness, temp, day, cosm)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(Richness, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "Richness",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "Richness", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = Richness, wid = cosm, 
  between = temp, within = day)
res.aov
# sphericity not violated, only main effect of day -> decrease over time


## Evenness
# Check assumptions
pan2 <- div %>% select(Evenness, temp, day, cosm)
pan2$day <- as.factor(pan2$day)

# Summary
pan2 %>%
  group_by(temp, day) %>%
  get_summary_stats(Evenness, type = "mean_sd")

bxp <- ggboxplot(
  pan2, x = "day", y = "Evenness",
  color = "temp", palette = temp_pal)
bxp

# Check normality
ggqqplot(pan2, "Evenness", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan2, dv = Evenness, wid = cosm, 
  between = temp, within = day)
res.aov
# sphericity not violated, effect of time and interaction


#### BARGRAPHS ####
### Bargraph on class level over time per treatment
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(ps, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class <- df2 %>% select(Sample, Class, Abundance, temp, day, rep)

# Rename classes for prettier plotting
class <- class %>%
  mutate(Class = recode(Class, "Prymnesiophyceae" = "Haptophyta")) %>%
  mutate(Class = recode(Class, "MOCH-1" = 'MOCH clade')) %>%
  mutate(Class = recode(Class, "MOCH-3" = 'MOCH clade')) %>%
  mutate(Class = recode(Class, "Haptophyta_Clade_HAP3" = 'Haptophyta')) %>%
  mutate(Class = recode(Class, "Prasino-Clade-VIII" = 'Prasinodermophyceae'))

## Classes of PR2 are spanning several taxonomic levels, rename to group
colnames(class)[2] <- "Group"

## Create color palette
class_pal <- qualpal(20, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting
class_plot <- ggplot(class, aes(fill = Group, x = day, y = Abundance)) +
  facet_wrap(rep ~ temp, ncol = 3) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_x_continuous(breaks = seq(0,27,3))+
  scale_fill_manual(values = class_pal$hex)

class_plot
ggsave("Output/PhytoGroups.png", class_plot, height = 10, width = 18, dpi = 320)

### Bargraph of dominant genera over time
genus <- phyloseq::tax_glom(ps, "Genus")
df <- plot_bar(genus, fill = "Genus")
df2 <- df$data
genus <- df2 %>% select(Sample, Genus, Abundance, temp, day, rep)

# Rename species for prettier plotting
genus$Genus[genus$Abundance < 100] <- "Other"
genus <- genus %>%
  mutate(Genus = recode(Genus, "Prymnesiophyceae_Clade_B4_X" = 'Prymnesiophyceae indet.')) %>%
  mutate(Genus = recode(Genus, "Parmales_env_1_X" = 'Parmales')) %>%
  mutate(Genus = recode(Genus, "Pedinellales_X" = 'Pedinellales'))

# Subsetting
genus_all <- genus
genus <- subset(genus, day > 12)

## Genera over whole time
# Create color palette
gen_pal_all <-
  qualpal(33, colorspace = list(
    h = c(0, 360),
    s = c(0.3, 1),
    l = c(0.2, 0.8)
  ))

leg_all <-
  c(
    "Aureococcus",
    "Bathycoccus",
    "Brockmanniella",
    "Chaetoceros",
    "Chrysochromulina",
    "Corethron",
    "Cylindrotheca",
    "Dictyocha",
    "Ditylum",
    "Florenciella",
    "Gephyrocapsa",
    "Haptolina",
    "Leptocylindrus",
    "Micromonas",
    "Minidiscus",
    "Odontella",
    "Ostreococcus",
    "Other",
    "Parmales",
    "Pedinellales",
    "Phaeocystis",
    "Picochlorum",
    "Plagioselmis",
    "Prymnesiophyceae indet.",
    "Prymnesium",
    "Pseudo-nitzschia",
    "Pterosperma",
    "Pyramimonas",
    "Skletonema",
    "Teleaulax",
    "Thalassionema",
    "Thalassiosira"
  )

## Plotting
genus_plot_all <- ggplot(genus_all, aes(fill = Genus, x = day, y = Abundance)) +
  facet_wrap(rep ~ temp, ncol = 3) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_x_continuous(breaks = seq(0,27,3)) +
  scale_fill_manual(labels = toexpr(leg_all,
                                    plain = c('Other', 'Prymnesiophyceae indet.')),
                    values = gen_pal_all$hex) +
  theme(legend.text.align = 0)

genus_plot_all
ggsave("Output/PhytoGeneraAll.png", genus_plot_all, height = 10, width = 16, dpi = 320)

## Genera only over experiment
# Create color palette
gen_pal <- qualpal(30, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

leg <- c("Bathycoccus", "Chaetoceros", "Chrysochromulina", "Corethron", "Cylindrotheca",
         "Dictyocha", "Ditylum", "Florenciella", "Gephyrocapsa", "Haptolina", "Leptocylindrus", "Micromonas",
         "Minidiscus", "Odontella", "Other", "Parmales", "Pedinellales", "Phaeocystis",
         "Picochlorum", "Plagioselmis", "Prymnesiophyceae indet.", "Prymnesium",
         "Pseudo-nitzschia", "Pterosperma", "Pyramimonas", "Skletonema", "Teleaulax", "Thalassionema",
        "Thalassiosira")

## Plotting
genus_plot <- ggplot(genus, aes(fill = Genus, x = day, y = Abundance)) +
  facet_wrap(rep ~ temp, ncol = 3) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_x_continuous(breaks = seq(0,27,3)) +
  scale_fill_manual(labels = toexpr(leg,
                                    plain = c('Other', 'Prymnesiophyceae indet.')),
                                    values = gen_pal$hex) +
  theme(legend.text.align = 0)

genus_plot
ggsave("Output/PhytoGeneraMain.png", genus_plot, height = 10, width = 16, dpi = 320)

### Bargraph of dominant species over time
species <- phyloseq::tax_glom(ps, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance, temp, day, rep)

# Rename species for prettier plotting
species$Species[species$Abundance < 100] <- "Other"
species <- species %>%
  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema')) %>%
  mutate(Species = recode(Species, "Chaetoceros_didymus_2" = 'Chaetoceros_didymus')) %>%
  mutate(Species = recode(Species, "Chaetoceros_debilis_1" = 'Chaetoceros_debilis')) %>%
  mutate(Species = recode(Species, "Chaetoceros_lorenzianus_1" = 'Chaetoceros_lorenzianus')) %>%
  mutate(Species = recode(Species, "Chaetoceros_lorenzianus_2" = 'Chaetoceros_lorenzianus')) %>%
  mutate(Species = recode(Species, "Micromonas_commoda_A2" = 'Micromonas_commoda')) %>%
  mutate(Species = recode(Species, "Prymnesiophyceae_Clade_B4_X_sp." = 'Prymnesiophyceae_indet.')) %>%
  mutate(Species = recode(Species, "Chrysophyceae_Clade-F_X_sp." = 'Chrysophyceae_indet.')) %>%
  mutate(Species = recode(Species, "Chrysophyceae_Clade-C_X_sp." = 'Chrysophyceae_indet.')) %>%
  mutate(Species = recode(Species, "Chrysophyceae_XXX_sp." = 'Chrysophyceae_indet.')) %>%
  mutate(Species = recode(Species, "Chaetoceros_sp_Clade_Na11C3" = 'Chaetoceros_sp.')) %>%
  mutate(Species = recode(Species, "Parmales_env_1_X_sp." = 'Parmales_indet.')) %>%
  mutate(Species = recode(Species, "Pedinellales_X_sp." = 'Pedinellales_indet.'))
species$Species <- gsub('_', ' ', species$Species)

## Create color palette
spe_pal <- qualpal(55, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting
species_plot <- ggplot(species, aes(fill = Species, x = day, y = Abundance)) +
  facet_wrap(rep ~ temp, ncol = 3) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  theme(axis.text.x = element_text(size = 14, face= "bold"),
        legend.text = element_text(face = "italic"))+
  scale_x_continuous(breaks = seq(0,27,3))+
  scale_fill_manual(values = spe_pal$hex)

species_plot
ggsave("Output/PhytoSpecies.png", species_plot, height = 10, width = 20, dpi = 320)

#### ORDINATION ####
### Compare different ordination methods

## Prepare datasets
ps_exp_raw <- subset_samples(ps_raw, day != "0" & day != "3" & day != "6" & day != "9" & day != "12")

## Ordinate
ordi_clr_euc <- ps_exp_raw %>%
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(shape = "day", size = 2) +
  scale_fill_manual(values=temp_pal) +
  scale_color_manual(values=uni_pal) +
  ggforce::geom_mark_ellipse(aes(fill = temp, color = uni_ID))

ordi_clr_euc

ordi_aitch <- ps_exp_raw %>%
  dist_calc("aitchison") %>%
  ord_calc() %>%
  ord_plot(shape = "day", size = 2) +
  scale_fill_manual(values=temp_pal) +
  scale_color_manual(values=uni_pal) +
  ggforce::geom_mark_ellipse(aes(fill = temp, color = uni_ID))

ordi_aitch

# Pattern: higher dispersion at 18°C than at the other two temperatures, at which it is comparable

## Save the simplest version
ggsave("Output/SpreadReplicates.png", ordi_aitch, dpi = 300, width = 8, height = 4)

#### BETADISPERSION ####
ps_var <- tax_glom(ps_clr, taxrank="Species")

# Subset for timepoints
ps_t15 <- subset_samples(ps_var, day =="15")
ps_t18 <- subset_samples(ps_var, day =="18")
ps_t21 <- subset_samples(ps_var, day =="21")
ps_t24 <- subset_samples(ps_var, day =="24")
ps_t27 <- subset_samples(ps_var, day =="27")

# Calculating aitchinson distances (euclidean of clr data) and beta-dispersion for each timepoint
aitch_15 <- distance(ps_t15, method = "euclidean")
samdf_15 <- data.frame(sample_data(ps_t15))
beta15 <- betadisper(aitch_15, samdf_15$temp)
dist15 <- data.frame(samdf_15)
dist15$betadisper <- beta15$distances

aitch_18 <- distance(ps_t18, method = "euclidean")
samdf_18 <- data.frame(sample_data(ps_t18))
beta18 <- betadisper(aitch_18, samdf_18$temp)
dist18 <- data.frame(samdf_18)
dist18$betadisper <- beta18$distances

aitch_21 <- distance(ps_t21, method = "euclidean")
samdf_21 <- data.frame(sample_data(ps_t21))
beta21 <- betadisper(aitch_21, samdf_21$temp)
dist21 <- data.frame(samdf_21)
dist21$betadisper <- beta21$distances

aitch_24 <- distance(ps_t24, method = "euclidean")
samdf_24 <- data.frame(sample_data(ps_t24))
beta24 <- betadisper(aitch_24, samdf_24$temp)
dist24 <- data.frame(samdf_24)
dist24$betadisper <- beta24$distances

aitch_27 <- distance(ps_t27, method = "euclidean")
samdf_27 <- data.frame(sample_data(ps_t27))
beta27 <- betadisper(aitch_27, samdf_27$temp)
dist27 <- data.frame(samdf_27)
dist27$betadisper <- beta27$distances

# create dataframe with distances of all timepoints
distances <- rbind(dist15, dist18, dist21, dist24, dist27)
rownames(distances) <- NULL
#distances <- na.omit(distances)
#distances$day <- as.factor(distances$day)
#distances <- as.data.frame(distances)

# Save as excel files
#write_xlsx(distances, "Data/BetadisperDistances_DNA.xlsx")

# Plot betadispersions over time per temperature
dist_plot <- ggplot(distances, aes(x=day, y=betadisper, color=temp, group=interaction(temp,day))) + 
  geom_boxplot(aes(color=temp, fill = temp), width = 2.25, size = 1)+
  geom_point(position = position_dodge(width=2.25))+
  labs(x="Incubation time (d)", y=bquote("Beta-dispersion")) +
  plot.theme+
  ggtitle("Species-level variation") +
  scale_x_continuous(breaks = seq(15, 27, 3))+
  scale_color_manual(values=temp_pal)+
  scale_fill_manual(values=temp_pal2)

dist_plot

#ggsave("Output/Betadispersion.png", dist_plot, dpi = 300, width = 8, height = 4)

### Repeated measures ANOVAs
## Subsets of 6 & 12 and 6 & 18
pan612 <- subset(distances, temp != "18")
pan618 <- subset(distances, temp != "12")

## For 6 to 12
# Check assumptions
pan612$day <- as.factor(pan612$day)

# Summary
pan612 %>%
  group_by(temp, day) %>%
  get_summary_stats(betadisper, type = "mean_sd")

bxp <- ggboxplot(
  pan612, x = "day", y = "betadisper",
  color = "temp", palette = temp_pal)
bxp

# exclude day 24 due to too few datapoints
pan612 <- subset(pan612, day != "24")
pan612$day <- factor(pan612$day)

# Check normality
ggqqplot(pan612, "betadisper", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan612, dv = betadisper, wid = cosm, 
  between = temp, within = day)

res.aov

res.aov_beta_612 <- get_anova_table(res.aov)
# nothing significant

## For 6 to 18
# Check assumptions
pan618$day <- as.factor(pan618$day)

# Summary
pan618 %>%
  group_by(temp, day) %>%
  get_summary_stats(betadisper, type = "mean_sd")

bxp <- ggboxplot(
  pan618, x = "day", y = "betadisper",
  color = "temp", palette = temp_pal)
bxp

# exclude day 24 due to too few datapoints
pan618 <- subset(pan618, day != "24")
pan618$day <- factor(pan618$day)

# Check normality
ggqqplot(pan618, "betadisper", ggtheme = theme_bw()) +
  facet_grid(day ~ temp, labeller = "label_both")
# very few datapoints, but apparently normal

# Compute 2-way RM ANOVA
res.aov <- anova_test(
  data = pan618, dv = betadisper, wid = cosm, 
  between = temp, within = day)

res.aov

res.aov_beta_618 <- get_anova_table(res.aov)
# significant effect of temperature and no interaction

#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
require(NCmisc)
packages <- list.functions.in.file("~/AWI/RProjects/TopTrons/CommunityComposition.R", alphabetic = TRUE) # set to your filepath
summary(packages)
