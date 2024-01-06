### TOPTRONS COMMUNITY COMPOSITION MICROZOOPLANKTON ###
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
setwd("~/AWI/RProjects/TopTrons/Supplement") # set to your working directory

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

# Create taxonomy file for microzooplankton
# based on https://doi.org/10.1111/jeu.12691
tax_chry <- filter(tax, tax$Class=="Chrysophyceae")
tax <- filter(tax, tax$Supergroup!="Archaeplastida")
tax <- filter(tax, tax$Division!="Haptophyta")
tax <- filter(tax, tax$Division!="Ochrophyta")
tax <- filter(tax, tax$Division!="Chlorophyta")
tax <- filter(tax, tax$Division!="Rhodophyta")
tax <- filter(tax, tax$Division!="Cryptophyta")
tax <- filter(tax, tax$Division!="Dinoflagellata")
tax <- filter(tax, tax$Division!="Prasinodermophyta")
tax <- filter(tax, tax$Division!="Chrompodellids")
tax <- filter(tax, tax$Class!="Chlorarachniophyceae") 

unique(tax$Division)
unique(tax$Class)

tax <- rbind(tax, tax_chry)

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
#rarecurve(t(otu_table(ps_raw)), step=50, cex=0.5) # looks fine
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

#### BARGRAPHS ####
### Bargraph on class level over time per treatment
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(ps, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class <- df2 %>% select(Sample, Class, Abundance, temp, day, rep)
class <- subset(class, day > 12)

# Rename classes for prettier plotting
class$Class[class$Abundance < 100] <- "Other"
class <- class %>%
  mutate(Class = recode(Class, "MAST-1" = "MAST clades")) %>%
  mutate(Class = recode(Class, "MAST-2" = "MAST clades")) %>%
  mutate(Class = recode(Class, "MAST-3" = "MAST clades")) %>%
  mutate(Class = recode(Class, "MAST-6" = "MAST clades")) %>%
  mutate(Class = recode(Class, "MAST-7" = "MAST clades")) %>%
  mutate(Class = recode(Class, "Apusomonadidae_Group-1" = 'Apusomonadidae')) %>%
  mutate(Class = recode(Class, "Centroheliozoa_X" = 'Centroheliozoa')) %>%
  mutate(Class = recode(Class, "Cercozoa_X" = 'Cercozoa')) %>%
  mutate(Class = recode(Class, "Picozoa_X" = 'Picozoa')) %>%
  mutate(Class = recode(Class, "Telonemia_X" = 'Telonemia')) %>%
  mutate(Class = recode(Class, "Filosa-Granofilosea" = 'Granofilosea')) %>%
  mutate(Class = recode(Class, "Filosa-Imbricatea" = 'Imbricatea')) %>%
  mutate(Class = recode(Class, "Filosa-Thecofilosea" = 'Thecofilosea')) %>%
  mutate(Class = recode(Class, "Pirsonia_Clade" = 'Pirsonia'))

## Classes of PR2 are spanning several taxonomic levels, rename to group
colnames(class)[2] <- "Group"

## Create color palette
class_pal <- qualpal(21, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting
class_plot <- ggplot(class, aes(fill = Group, x = day, y = Abundance)) +
  facet_wrap(rep ~ temp, ncol = 3) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_x_continuous(breaks = seq(0,27,3))+
  scale_fill_manual(values = class_pal$hex)

class_plot
ggsave("Output/MicrozooGroups.png", class_plot, height = 10, width = 18, dpi = 320)
