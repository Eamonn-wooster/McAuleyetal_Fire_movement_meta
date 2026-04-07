## Adding Body Mass Data from ReptTraits (https://figshare.com/articles/dataset/ReptTraits_a_comprehensive_database_of_ecological_traits_in_reptiles/24572683) 
## and AvoNet (https://figshare.com/s/b990722d72a26b5bfead)
## and Phylacine (https://github.com/MegaPast2Future/PHYLACINE_1.2/tree/v1.2.1)
## CM

## Packages ##

library(stringr)
library(tidyverse)
library(here)
library(dplyr)

#Setting Directory

here()


################################# Reptiles ###################################

#Reading in extraction data
dat <- read.csv(here("Outputs/Dat_for_analysis-ss-tree.csv")) 

# Reading in ReptTraits v1.2
repttraits <- read.csv(here("Data/Traits/ReptTraits_dat.csv"))

dat <- dat %>%
  rename(Species = Study_Species)

combined_data_rept <- dat %>%
  left_join(repttraits %>% select(Species, Maximum_body_mass_g), by = "Species")

# Renaming column to include source and descriptor

combined_data_rept = combined_data_rept %>% 
  rename(ReptTraits.Max_body_mass_g = Maximum_body_mass_g)

# adding value for Terrapene carolina carolina based on Terrapene carolina mass

combined_data_rept$ReptTraits.Max_body_mass_g[dat$Species == "Terrapene carolina carolina"] <- 1093.0

####################################### Aves ##################################
# Reading in AvoNet

avonet = read.csv(here("Data/Traits/AVONET1_BirdLife.csv"))

# Renaming Species1 to match column name in dat

avonet <- avonet %>%
  rename(Species = Species1)

combined_data_rept.aves <- combined_data_rept %>%
  left_join(avonet %>% select(Species, Mass), by = "Species")

## renaming aves mass to include source and descriptor

combined_data_rept.aves <- combined_data_rept.aves %>%
  rename(AvoNet.Mean_mass_g = Mass)

## Mutating the ReptTraits and Avonet body mass values into a single column

combined_data_rept.aves <- combined_data_rept.aves %>%
  mutate(BodyMass.g = coalesce(ReptTraits.Max_body_mass_g, AvoNet.Mean_mass_g))

## Saving a new csv file with body mass 

 # write_csv(combined_data_rept.aves, "Outputs/Dat_for_analysis-ss-aves-rep-mass-Rev3.csv")

######################## Mammals ##############################################

# reading in extraction data, with ReptTraits and Avonet already added

 # dat <- read.csv(here("Outputs/Dat_for_analysis-ss-aves-rep-mass-Rev3.csv"))

# Reading in Phylacine data

phy = read.csv(here("Data/Traits/Phylacine_Trait_data.csv"))

# Renaming Binomial.1.2 to match column name in dat

phy <- phy %>%
  rename(Species_tree = Binomial.1.2)

# adding body mass data

combined_data_rept.aves.mam <- combined_data_rept.aves %>%
  left_join(phy %>% select(Species_tree, Mass.g), by = "Species_tree")

# Renaming column to include source and descriptor

combined_data_rept.aves.mam = combined_data_rept.aves.mam %>% 
  rename(Phy.body_mass_g = Mass.g)

# Coalescing body mass data into a single column 

combined_data_rept.aves.mam <- combined_data_rept.aves.mam %>%
  mutate(BodyMass.g = coalesce(ReptTraits.Max_body_mass_g, AvoNet.Mean_mass_g, Phy.body_mass_g))

# Checking for NA's, to add manually

summary(combined_data_rept.aves.mam$BodyMass.g)

table(is.na(combined_data_rept.aves.mam$BodyMass.g))

# Writing csv

write_csv(combined_data_rept.aves.mam, "Outputs/Dat_for_analysis-ss-aves-rep-mam-mass-Rev4.csv")
  
# Dat_for_analysis-ss-aves-rep-mam-mass-Rev4.csv updated to Dat_for_analysis-mass-Rev4 after manual updates 

