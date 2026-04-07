####################### Adding Order and Class to Data #######################

# Packages

library(taxize)

pacman::p_load(devtools, 
               tidyverse, 
               R.rsp, 
               dplyr,
               broom,
               tidyverse,
               here
               )

######## Extracting class ####################################################

# Setting Directory

here()

#Importing data

dat <- read.csv(here("Outputs/Dat_for_analysis-ss-mass-biome-clim-fireact-Rev4.csv"))

# Extracting Species

species_list = dat$Species_tree

# Getting Class and Order

class_df <- tax_name(species_list, get = "class", db = "ncbi")
order_df <- tax_name(species_list, get = "order", db = "ncbi")

# Renaming for Join

order_df = rename(order_df, Species_tree = query)
class_df = rename(class_df, Species_tree = query)

# Making sure class_df has only one row per species, to avoid unwanted duplication

class_df_unique <- class_df %>%
  distinct(Species_tree, class)

# Join

dat <- dat %>%
  left_join(class_df_unique, by = "Species_tree")

# Making sure order_df has only one row per species, to avoid unwanted duplication

order_df_unique <- order_df %>%
  distinct(Species_tree, order)

# Join

dat <- dat %>%
  left_join(order_df_unique, by = "Species_tree")

# Naming columns to include database source

dat = rename(dat, class_ncbi = class)
dat = rename(dat, order_ncbi = order)

# Writing CSV

write.csv(dat, here("Outputs/Dat_for_analysis-ss-mass-biome-clim-fireact-tax-Rev4.csv"))
