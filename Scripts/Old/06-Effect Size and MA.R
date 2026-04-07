
# install.packages("pacman")

#devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
# pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans,
#               ape, phytools, flextable)

rm(list = ls())

### Loading packages ----
pacman::p_load(devtools, 
               tidyverse, 
               metafor, 
               patchwork, 
               R.rsp, orchaRd,
               emmeans,
               metafor,
               orchaRd,
               stringr,
               dplyr,
               broom,
               tidyverse,
               here,
               broom,
               ggplot2,
               data.table,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp #,#added
               # miWQS #added
) 

## Setting WD

here()

dat <- read.csv(here("Outputs/Dat_for_analysis-ss-mass-biome-tax-Rev3.csv")) 

# getting tree - making new trees for each evidence catergories has no influence on the model

tree1 <- read.tree(here("Outputs/Tree_fire-Rev3.tre")) 

# getting branch length and correlation matrix

tree1b <- compute.brlen(tree1)
cor1 <-  vcv(tree1b, corr=T)

#' [EW - this column is in some dataframes and not others - make sure when you add new data/change the data set you carry this across]
#' [ i have been manually adding this bc the tree and dataset are struggling to talk to each other , not sure why]

## checking the match 

setdiff(dat$Species_tree, tree1$tip.label)

#' [EW - This HAS to return character(0) for the model to run]
#' [ EW if it does not you haven't carried over the species tree column - go back and do this]


# # creating non-phylo columns - has to be tree column 
dat$Species2 <- dat$Species_tree

#################################################################################
######################### Converting Error ######################################
#################################################################################

### >>> functions for converting error ------

#Function for converting CI's to SE
ci_to_se <-function(upper, lower, n){
  # upper: upper limit of CI
  # lower: lower limit of CI
  # n: sample size
  # return: SE
  se_value<-(upper - lower) / (2 * qt(0.975, n-1))
  return(se_value)
}

#function for converting SE to SD
se_to_sd <- function(se, n) {
  # se: standard error
  # n: sample size
  # return: standard deviation
  sd_value <- se * sqrt(n)
  return(sd_value)
}

#converting error types

unique(dat$CI_Type)

### >>> Standard error --------

SD <- filter(dat, CI_Type == "SD")

#Taking the mean from the upper gives us the SD values given that SD is a single number (same applies to SE below)

SD$sd_Unburned <- SD$Upper_Unburned - SD$Mean_Unburned

SD$sd_Fire <- SD$Upper_Fire - SD$Mean_Fire

### >>> 95% CI's --------

CI <- filter(dat, CI_Type == "CI")

#Squaring gets it from SE to SD

CI$sd_Unburned <- ci_to_se(CI$Upper_Unburned, CI$Lower_Unburned, CI$Unburned_n)^2

CI$sd_Fire <- ci_to_se(CI$Upper_Fire, CI$Lower_Fire, CI$Burned_n)^2

### >>> Standard deviation --------

SE <- filter(dat, CI_Type == "SE")

# Taking the upper value from the mean gives us the SD values

SE$sd_Unburned <- SE$Upper_Unburned - SE$Mean_Unburned

SE$sd_Fire <- SE$Upper_Fire - SE$Mean_Fire

# now convert SE to SD

SE$sd_Unburned <- se_to_sd(SE$sd_Unburned, SE$Unburned_n)

SE$sd_Fire <- se_to_sd(SE$sd_Fire, SE$Burned_n)

# Put them back together!! 

data <- rbind(SD, CI, SE)

#################################################################################
############################ Cleaning Data ######################################
#################################################################################

## Tidying up movement category 

# looking at unique category counts to identify duplicates 

table(data$Movement_category)

# Standardising all to lower case

data$Movement_category <- tolower(trimws(data$Movement_category))

# Back to Title case

data$Movement_category <- str_to_title(data$Movement_category)

# To factor, to use in analysis

data$Movement_category <- factor(data$Movement_category)

table(data$Movement_category) ## good to go :)

## Cleaning Fire Category

# looking at unique category counts to identify duplicates 

table(data$Impact_fire_type)

data$Impact_fire_type <- factor(data$Impact_fire_type)

data$class_ncbi = factor(data$class_ncbi)
table(data$class_ncbi)

#################################################################################
############################ Overall Effect Size ################################
#################################################################################

# Calculate Hedges' g

es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = data,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(es) ## CM yi is the SMD effect size, and vi is the sampling variance. 

#The variance here seems better?? I think so too - EW

#obs level random effect

es$Obs_ID <- factor(1:nrow(es))

#This was something strange replacing with correct

es$StudyID <- NULL

es$Study_ID <- factor(as.numeric(as.factor(es$Firstauthor_Year)))


#################################################################################
################################## MA Models ####################################

### Overall models

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                       data = es) 


if(rerun){
  mod.overall <- rma.mv(yi = yi, V = vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  es,
                   #mods = XXX - uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall, here("Outputs", "mod.overall.Rds"))
  summary(mod.overall)
}else{
  mod.overall <- readRDS(here("Outputs", "mod.overall.Rds"))
}

overall <- orchard_plot(mod.overall, xlab = "Change in movement (Hedge's g)", group = "Study_ID")
overall


###Movement categories 

#' [EW - we don't need to run movement categories as seperate models - you can but this is also ok. I see i told you to do that - sorry!]
#' [EW - honestly it works either way, but maybe we want to compare between types of movement, if so, run a single model]
#' [sorry, bad advice on running seperate models, didn't fully understand the data then :/ ]
if(rerun){
  mod.overall_move <- rma.mv(yi = yi, V = vcv,
                        random = list(~1 | Study_ID,
                                      ~1 | Species_tree, # phylo effect 
                                      ~1 | Species2, # non-phylo effect 
                                      ~1 | Obs_ID), 
                        data =  es,
                        mods = ~ Movement_category -1,
                        control=list(optimizer="BFGS"),
                        test = "t",
                        sparse = TRUE,
                        R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall_move, here("Outputs", "mod.overall_move.Rds"))
  summary(mod.overall_move)
}else{
  mod.overall_move <- readRDS(here("Outputs", "mod.mod.overall_move.Rds"))
}

overall_move <- orchard_plot(mod.overall_move, xlab = "Change in movement (Hedge's g)", group = "Study_ID", mod = "Movement_category",
                             angle = 45)
overall_move

#################################################################################
##################### Movement Category Models: Home Range ######################
#################################################################################

# Filtering Home Range

homerange.dat <- data %>%
  filter(Movement_category == "Home Range")

# calculating effect size

homerange.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = homerange.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(homerange.es)

homerange.es$Obs_ID <- factor(1:nrow(homerange.es))

#This was something strange replacing with correct

homerange.es$StudyID <- NULL

homerange.es$Study_ID <- factor(as.numeric(as.factor(homerange.es$Firstauthor_Year)))

# VCV

homerange.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = homerange.es) 

# Overall Home Range Model

#' [EW - this should be false so the model doesn't run each time]
rerun <- FALSE

if(rerun){
  mod.hr <- rma.mv(yi = yi, V = homerange.vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species_tree, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  homerange.es,
                 #mods = XXX - uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                 test = "t",
                 sparse = TRUE,
                 R = list(Species_tree = cor1)
  )
  saveRDS(mod.hr, here("Outputs", "mod_hr.Rds"))
  summary(mod.hr)
}else{
  mod.hr <- readRDS(here("Outputs", "mod_hr.Rds"))
}

# Plotting and viewing model results

homerange <- orchard_plot(mod.hr, xlab = "Change in Home Range Size (Hedge's g)", group = "Study_ID")
homerange

#Saving plot

ggsave(here("Outputs", "Figures", "homerange_test.pdf"), width = 6, height = 4)

##################################################################################
################# Movement Category Models: Linear Disatnce ######################
##################################################################################

# Filtering Distance Data

table(data$Movement_category)

ld.dat <- data %>%
  filter(Movement_category %in% c("Linear Distance", "Linear Distance - Speed"))

# calculating effect size

ld.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = ld.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(ld.es)

ld.es$Obs_ID <- factor(1:nrow(ld.es))

#This was something strange replacing with correct

ld.es$StudyID <- NULL

ld.es$Study_ID <- factor(as.numeric(as.factor(ld.es$Firstauthor_Year)))

# VCV

ld.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                       data = ld.es) 

# Overall Home Range Model


if(rerun){
  mod.ld <- rma.mv(yi = yi, V = ld.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  ld.es,
                   #mods = XXX - uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.ld, here("Outputs", "mod_ld.Rds"))
  summary(mod.ld)
}else{
  mod.ld <- readRDS(here("Outputs", "mod_ld.Rds"))
}

# Plotting and viewing model results

lineardist <- orchard_plot(mod.ld, xlab = "Change in Linear Distance (Hedge's g)", group = "Study_ID")
lineardist

#Saving plot

ggsave(here("Outputs", "Figures", "lineardist_test.pdf"), width = 6, height = 4)

## CM: Need to check this model against a model without Linear Distance - Speed
# to ensure validity of grouping 

#################################################################################
################# Movement Category Models: Space Use ###########################
#################################################################################

# Filtering Space use data

table(data$Movement_category)

su.dat <- data %>%
  filter(Movement_category %in% c("Other Space Use - Proportion", "Other Space Use - Time"))

# calculating effect size

su.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = su.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(su.es)

su.es$Obs_ID <- factor(1:nrow(su.es))

#This was something strange replacing with correct

su.es$StudyID <- NULL

su.es$Study_ID <- factor(as.numeric(as.factor(su.es$Firstauthor_Year)))

# VCV

su.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                data = su.es) 

# Overall Home Range Model


if(rerun){
  mod.su <- rma.mv(yi = yi, V = su.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  su.es,
                   #mods = XXX - uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.su, here("Outputs", "mod_su.Rds"))
  summary(mod.su)
}else{
  mod.su <- readRDS(here("Outputs", "mod_su.Rds"))
}

# Plotting and viewing model results

spaceuse <- orchard_plot(mod.su, xlab = "Change in Space Use (Hedge's g)", group = "Study_ID")
spaceuse

#Saving plot

ggsave(here("Outputs", "Figures", "spaceuse)test.pdf"), width = 6, height = 4)

## CM: Need to check this model against a model without one of the space use groups
# not sure how this would go due to small sample sizes

#################################################################################
################ Movement Category Models: Tortuosity  ##########################
#################################################################################

# Filtering Tort data

table(data$Movement_category)

tort.dat <- data %>%
  filter(Movement_category %in% c("Tortuosity"))

# calculating effect size

tort.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = tort.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(tort.es)

tort.es$Obs_ID <- factor(1:nrow(tort.es))

#This was something strange replacing with correct

tort.es$StudyID <- NULL

tort.es$Study_ID <- factor(as.numeric(as.factor(tort.es$Firstauthor_Year)))

# VCV

tort.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                data = tort.es) 

# Overall Home Range Model

if(rerun){
  mod.tort <- rma.mv(yi = yi, V = tort.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   control=list(optimizer="BFGS"), #This only if model doesn't run
                   data =  tort.es,
                   #mods = XXX - uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.tort, here("Outputs", "mod_tort.Rds"))
  summary(mod.tort)
}else{
  mod.tort <- readRDS(here("Outputs", "mod_tort.Rds"))
}

# Plotting and viewing model results

tort <- orchard_plot(mod.tort, xlab = "Change in Tortuosity (Hedge's g)", group = "Study_ID")
tort

#Saving plot

ggsave(here("Outputs", "Figures", "tort_test.pdf"), width = 6, height = 4)

# CM: Probably insufficent data to model this
#' [EW - yeah there isn't much but thats ok]



###############################################################################
############################# Moderator: Fire Type ############################
###############################################################################

# Filtering out not provided, unclear or both fire types

#' [EW - this is fine but you don't need to recalc the effect size each time, instead you can filter the es dataset]
table(data$Impact_fire_type)

fire.dat <- data %>%
  filter(Impact_fire_type %in% c("Planned", "Wildfire"))

fire.dat$Impact_fire_type <- factor(fire.dat$Impact_fire_type)

# calculating effect size

fire.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = fire.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(fire.es)

fire.es$Obs_ID <- factor(1:nrow(fire.es))

#This was something strange replacing with correct

fire.es$StudyID <- NULL

fire.es$Study_ID <- factor(as.numeric(as.factor(fire.es$Firstauthor_Year)))

# VCV

fire.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                       data = fire.es) 
# Trimming Tree

#' [EW - neat idea, i note we get a warning but if you compare the models, having empty values in the correlation just drops them, making the model identical]
#' [EW - This means you don't need to trim the tree bc we get the same results either way - please let me know if you find it changes your results]
#' [EW - i reran you model with cor1 and cor1_fire and they were identical]
#' [EW - i'm not sure what the redundent predictor error is but the model seems ok as does the set up]
species_used <- unique(fire.es$Species_tree)

cor1_fire <- cor1[species_used, species_used]

# Model with Fire Type


if(rerun){
  mod.fire <- rma.mv(yi = yi, V = fire.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  fire.es,
                   mods = ~ Impact_fire_type,
                   control=list(optimizer="BFGS"),
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1_fire)
  )
  saveRDS(mod.fire, here("Outputs", "mod_fire.Rds"))
  summary(mod.fire)
}else{
  mod.fire <- readRDS(here("Outputs", "mod_fire.Rds"))
}

plot.fire <- orchard_plot(mod.fire, xlab = "Change in movement (Hedge's g)", group = "Impact_fire_type", mod = "Impact_fire_type")
plot.fire

ggsave(here("Outputs", "Figures", "Fire_type.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 
round(i2_ml(mod.fire), 2) #' [You only need this from the overall model - we just report I2 onces]

###############################################################################
######################## Model: Class as Moderator ############################
###############################################################################

table(data$class_ncbi)

# Filtering out Amphibia due to small sample size - agreed! 

class.dat <- data %>%
  filter(class_ncbi %in% c("Aves", "Mammalia", "Reptilia"))

class.dat$class_ncbi <- factor(class.dat$class_ncbi)

# Effect Size

class.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = class.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(class.es)

class.es$Obs_ID <- factor(1:nrow(class.es))

#This was something strange replacing with correct

class.es$StudyID <- NULL

class.es$Study_ID <- factor(as.numeric(as.factor(class.es$Firstauthor_Year)))

# VCV

class.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                  data = class.es) 

# Trimming tree matrix 
#' [EW - same here no difference in model]
species.class <- unique(class.es$Species_tree)

cor1_class <- cor1[species.class, species.class]

# Model
if(rerun){
  mod.class <- rma.mv(yi = yi, V = class.vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species_tree, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  class.es,
                 mods = ~ class_ncbi -1, #remove intercept because we want a value for each group, not comparison between each.
                   test = "t",
                 sparse = TRUE,
                 R = list(Species_tree = cor1_class)
  )
  saveRDS(mod.class, here("Outputs", "mod_class.Rds"))
  summary(mod.class)
}else{
  mod.class <- readRDS(here("Outputs", "mod_class.Rds"))
}

plot.class <- orchard_plot(mod.class, xlab = "Change in movement (Hedge's g)", group = "class_ncbi", mod = "class_ncbi")
plot.class

ggsave(here("Outputs", "Figures", "class.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 
round(i2_ml(mod.class), 2)

###############################################################################
######################### Biome as Moderator ##################################
###############################################################################

table(data$biome_name)

# Filtering out Montane grasslands and shrublands and Tropical and subtropical 
#coniferous forests due to small sample size

biome.dat <- data %>%
  filter(biome_name %in% c("Boreal forests/taiga ", 
                           "Deserts and xeric shrublands",
                           "Flooded grasslands and savannas", 
                           "Mediterranean forests, woodlands, and scrub or sclerophyll forests", 
                           "Temperate broadleaf and mixed forests",
                           "Temperate coniferous forests",
                           "Temperate grasslands, savannas, and shrublands",
                           "Tropical and subtropical grasslands, savannas, and shrublands",
                           "Tropical and subtropical moist broadleaf forests"))

biome.dat$biome_name <- factor(biome.dat$biome_name)

# Effect Size

biome.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = biome.dat,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(biome.es)

biome.es$Obs_ID <- factor(1:nrow(biome.es))

#This was something strange replacing with correct

biome.es$StudyID <- NULL

biome.es$Study_ID <- factor(as.numeric(as.factor(biome.es$Firstauthor_Year)))

# VCV

biome.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = biome.es) 

# Trimming tree matrix

species.biome <- unique(biome.es$Species_tree)

cor1_biome <- cor1[species.biome, species.biome]

# Model

if(rerun){
  mod.biome <- rma.mv(yi = yi, V = biome.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  biome.es,
                      mods = ~ biome_name -1, #' [We remove the intercept here bc we want estimates for all values and not a comparision to the first one]
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1_biome) #' [EW - this is the same as above]
  )
  saveRDS(mod.biome, here("Outputs", "mod_biome.Rds"))
  summary(mod.biome)
}else{
  mod.biome <- readRDS(here("Outputs", "mod_biome.Rds"))
}

plot.biome <- orchard_plot(mod.biome, xlab = "Change in movement (Hedge's g)", group = "biome_name", mod = "biome_name")

#'[EW - check this plot labels and work correctly since i removed the intercept from the model!]


#Fixing labels

wrapped_labels <- function(x) str_wrap(x, width = 25)

plot.biome + 
  scale_x_discrete(labels = wrapped_labels) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 8)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 8)) +
  theme(axis.title.x = element_text(size = 10))

ggsave(here("Outputs", "Figures", "Biome.pdf"), width = 8, height = 6)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 
round(i2_ml(mod.biome), 2)

###############################################################################
############################ Body Mass model ##################################
###############################################################################

#We need to log trans form BM before modelling

es$logmass <- log(es$BodyMass.g)

if(rerun){
  mod.bm <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species_tree, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ logmass, #uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                 test = "t",
                 sparse = TRUE,
                 R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm, here("Outputs", "mod1_bm.Rds"))
  summary(mod.bm)
}else{
  mod.bm <- readRDS(here("Outputs", "mod1_bm.Rds"))
}

mod.bm


#Extracting and plotting body mass models 

round(r2_ml(mod.bm)*100, 2) #This is r2 of the line - CM
#I've replaced mod1_func with mod1, is this correct? 
#' [EW - yes this is correct]


## replace logpredmass with you body mass column and herb_func with you data frame
newdat <- data.frame(logmass = seq(min(es$logmass, na.rm = TRUE),
                                       max(es$logmass, na.rm = TRUE),
                                       length.out = 100)) ## CM: Why is the length out 100?

X <- model.matrix(~ logmass, data = newdat)[, -1, drop = FALSE]

preds <- predict(mod.bm, newmods = X)

plot_data <- cbind(newdat, preds)

setDT(es)

es[, wi := 1/sqrt(vi)]
es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
es

prey_bm_text <- tidy(mod.bm)
prey_bm_text ## Very weak association of larger body mass = slightly smaller effect sizes

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.mean

ggsave(here("Outputs/Figures", "body mass.pdf"), height = 4, width = 8)


###############################################################################
############################ Treatment Age Model ##############################
###############################################################################

# Coercing variable to integer

data <- data %>%
  mutate(Maximum_treatment_burn_age_days = if_else(
    grepl("^[-]?[0-9]+$", Maximum_treatment_burn_age_days),
    as.integer(Maximum_treatment_burn_age_days),
    NA_integer_
  ))

# Looking for NAs

summary(data$Maximum_treatment_burn_age_days)

# Removing NAs

data.age <- data[!is.na(data$Maximum_treatment_burn_age_days), ]

# Effect Size

age.es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = data.age,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(age.es)

age.es$Obs_ID <- factor(1:nrow(age.es))

#This was something strange replacing with correct

age.es$StudyID <- NULL

age.es$Study_ID <- factor(as.numeric(as.factor(age.es$Firstauthor_Year)))

# VCV

age.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = age.es) 

# Trimming Tree

species.age <- unique(age.es$Species_tree)

cor1_age <- cor1[species.age, species.age]

# Model

hist(age.es$Maximum_treatment_burn_age_days)

#' [EW - i'd try rerunning this after dropping the crazy outlier and then consider log transforming if you need to make it a but more normal]

if(rerun){
  mod.age <- rma.mv(yi = yi, V = age.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  age.es,
                   mods = ~ Maximum_treatment_burn_age_days, #uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1_age)
  )
  saveRDS(mod.age, here("Outputs", "mod1_age.Rds"))
  summary(mod.age)
}else{
  mod.age <- readRDS(here("Outputs", "mod1_age.Rds"))
}

mod.age


#Extracting and plotting body mass models 

round(r2_ml(mod.age)*100, 2) #This is r2 of the line - CM
#I've replaced mod1_func with mod1, is this correct? 


## ### replace logpredmass with you body mass column and herb_func with you data frame

newdat.age <- data.frame(Maximum_treatment_burn_age_days = seq(min(age.es$Maximum_treatment_burn_age_days, na.rm = TRUE),
                                      max(age.es$Maximum_treatment_burn_age_days, na.rm = TRUE),
                                      length.out = 100)) ## CM: Why is the length out 100?

X.age <- model.matrix(~ Maximum_treatment_burn_age_days, data = newdat.age)[, -1, drop = FALSE]

preds.age <- predict(mod.age, newmods = X.age)

plot_data.age <- cbind(newdat.age, preds.age)

setDT(age.es)

age.es[, wi := 1/sqrt(vi)]
age.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
age.es

prey_age_text <- tidy(mod.age)
prey_age_text

#' [EW - this should work better now, effect size should be on the y not se - lemme know if it doesn't plot like the bm one above]

p.age.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = plot_data.age, 
              aes(x = Maximum_treatment_burn_age_days, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.age,
              aes(x = Maximum_treatment_burn_age_days,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.age, 
              aes(x = Maximum_treatment_burn_age_days,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.age, 
            aes(x = Maximum_treatment_burn_age_days,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Maximum Burn Age (days)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.age.mean

ggsave(here("Outputs/Figures", "Burn Age.pdf"), height = 4, width = 8)


################################################################################
################################################################################

#------------------------------------------------------------------------------#
# >>> Overall model from EW ----------------------------------------------

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 


# for the model loops
rerun <- F

# >>> Overall model ----------------------------------------------
rerun <- TRUE

if(rerun){
  mod1 <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species_tree, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ##~ class_ncbi,
                   test = "t",
                 sparse = TRUE,
                 R = list(Species_tree = cor1)
  )
  saveRDS(mod1, here("Outputs", "mod1_overall.Rds"))
  summary(mod1)
}else{
  mod1 <- readRDS(here("Outputs", "mod1_overall.Rds"))
}

overall <- orchard_plot(mod1, xlab = "Change in movement (Hedge's g)", group = "Study_ID")

overall

ggsave(here("Outputs", "Figures", "overlap_overall_test.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 
round(i2_ml(mod1), 2)


# [EW - i'd recommend running each movement category as a separate model, i.e., 
# filtering for home range, making a new vcv, see above, just replace data with 
# new data in that code, and then running with no moderators, see what that gives 
# you and it allows you to run each category with body mass]
# [Also, lets clean these movement categories into something more manageable,
#  i.e., home range, linear distance, space use, tortuosity should be fine - this is a bit of a mess at the moment. ]

if(rerun){
  mod2 <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species_tree, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Movement_category,
                 control=list(optimizer="BFGS"), #Model was struggling with convergence
                 test = "t",
                 sparse = TRUE,
                 R = list(Species_tree = cor1)
  )
  saveRDS(mod2, here("Outputs", "mod2_move.Rds"))
  summary(mod2)
}else{
  mod2 <- readRDS(here("Outputs", "mod1_overall.Rds"))
}

move_cat <- orchard_plot(mod2, xlab = "Change in movement (Hedge's g)", group = "Study_ID")

overall


#' [ EW leaving you to do the rest - below is code for nice models from body mass models instead or Orchidplots, 
#' [they do continuous predictors real ugly]
#' [ i've also let the mammal body mass for you to do, since you are super capable, use Phylacine or Mass of Mammals]
#' [Make sure you don't lose any of the columns in this dataset]
#' 



####CM models below -------------------
## Running model ##

## Running base model, with defualt rma() settings ##

res <- rma(yi = yi, vi = vi, data = es, method = "REML")

# View summary results
summary(res)

#forest plot
forest(res, xlab = "Hedges' g")

#funnel plot
funnel(res)
regtest(res)

## body mass

res_mod <- rma(yi, vi, mods = ~ BodyMass.g, data = es, method = "REML")

# View results
summary(res_mod)


# Making sure movement category is a factor for use in analysis

es$Movement_category = factor(es$Movement_category)
str(es$Movement_category)

## Movement category ##

levels <- levels(es$Movement_category)
results <- lapply(levels, function(cat) {
  subset_data <- subset(es, Movement_category == cat)
  rma(yi, vi, data = subset_data, method = "REML")
})


# Print summaries
for (i in seq_along(levels)) {
  cat("\n\nSubgroup:", levels[i], "\n")
  print(summary(results[[i]]))
}

for (movement_cat in levels) {
  subset_data <- subset(es, Movement_category == movement_cat)
  
  # Check if there are enough studies in the subgroup
  if (nrow(subset_data) > 1) {
    res <- rma(yi, vi, data = subset_data, method = "REML")
    
    # Generate forest plot
    forest(res, main = paste("Forest Plot for", movement_cat), xlab = "Effect Size (Hedges' g)")
  } else {
    message("Subgroup '", movement_cat, "' has too few studies. Skipping.")
  }
}





