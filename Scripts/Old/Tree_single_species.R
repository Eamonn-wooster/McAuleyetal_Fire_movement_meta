#title: Phylo Tree for multi-level meta-analysis

#Author: EIF Wooster 


rm(list = ls())
library("pacman")
#install.packages("pacman")

#install.packages("remotes")
#remotes::install_github("ropensci/rotl", force = TRUE)

#install.packages("rotl")# You'll have to install the package below
pacman::p_load(devtools, 
               tidyverse, 
               R.rsp, 
               dplyr,
               broom,
               tidyverse,
               here,
               broom,
               data.table,
               rotl,
               ape)

#getwd()

#Read in your data
dat <- read.csv(here("data","revision", "Single_species_overlap_revision.csv")) 

##############################
#########Species tree############
##############################

#Pull out all species - leaving object as Species 
Species <- dplyr::select(dat, Species)#after column is you latin species name column

Species <- unique(Species)

length(Species$Species)

setDT(Species)

#This checks the species of yours that are not in the tree
is_in_tree(ott_id(Species$Species))

##Evything needs to be at the species level for the most part
# You can use the code examples below to change names of species

Species <- Species[Species == "Canis_lupus_arabs", Species := "Canis_lupus"]

Species <- unique(Species)

#Reaching out to tree to see matches - make sure all your species have a correct match
resolved_Species <- tnrs_match_names(Species$Species, context_name = "Animals")

resolved_Species

#length(resolved_Species$ott_id)

Species_tree <- tol_induced_subtree(ott_ids = resolved_Species$ott_id)

length(Species_tree$tip.label) 

# getting ride of the ott id   
Species_tree$tip.label <- gsub("_ott\\d+", "", Species_tree$tip.label)


setdiff(names_rtol, Species_tree$tip.label)
setdiff(Species_tree$tip.label, Species$Species)

tree1 <- Species_tree

plot(tree1, no.margin = TRUE)

# clearning up a bit


# writing tree
write.tree(tree1, file = here("Tree","Tree_fire.tre"))

