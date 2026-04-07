################### Effect Size & Meta analysis Modelling: No space Use #####################
## Authors: EW & CM


#################################################################################
######################### Set Up & Reading in Data ##############################
#################################################################################

# Loading packages 

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
               viridis,
               data.table,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp 
               #added
               # miWQS #added
) 

# Setting WD

here()

# Reading in Data

dat <- read.csv(here("Outputs/No_space_use/Dat_for_analysis-ss-mass-biome-clim-fireact-tax-Rev5.csv")) 

# Reading in tree

tree1 <- read.tree(here("Outputs/Tree_fire-Rev4.tre")) 

# Getting branch length and correlation matrix

tree1b <- compute.brlen(tree1)
cor1 <-  vcv(tree1b, corr=T)

# checking the match 

setdiff(dat$Species_tree, tree1$tip.label)

#' [EW - This HAS to return character(0) for the model to run]

# creating non-phylo columns - has to be tree column 

dat$Species2 <- dat$Species_tree

#################################################################################
######################### Converting Error ######################################
#################################################################################

### functions for converting error ------

# Function for converting CI's to SE

ci_to_se <-function(upper, lower, n){
  # upper: upper limit of CI
  # lower: lower limit of CI
  # n: sample size
  # return: SE
  se_value<-(upper - lower) / (2 * qt(0.975, n-1))
  return(se_value)
}

# Function for converting SE to SD

se_to_sd <- function(se, n) {
  # se: standard error
  # n: sample size
  # return: standard deviation
  sd_value <- se * sqrt(n)
  return(sd_value)
}

### Converting error types

unique(dat$CI_Type)

### >>> Standard deviation -------- 
#' [EW - the title for this section is confusing - i've renamed make sure they make sense, i think SD and SE were backwards]

SD <- filter(dat, CI_Type == "SD")

#Taking the mean from the upper gives us the SD values given that SD is a single number (same applies to SE below)

SD$sd_Unburned <- SD$Upper_Unburned - SD$Mean_Unburned

SD$sd_Fire <- SD$Upper_Fire - SD$Mean_Fire

### >>> 95% CI's --------

CI <- filter(dat, CI_Type == "CI")

#Squaring gets it from SE to SD

CI$sd_Unburned <- ci_to_se(CI$Upper_Unburned, CI$Lower_Unburned, CI$Unburned_n)^2

CI$sd_Fire <- ci_to_se(CI$Upper_Fire, CI$Lower_Fire, CI$Burned_n)^2

### >>> Standard error -------- 

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

# Cleaning Categories

table(data$Movement_category) 

data <- data %>%
  mutate(Movement_category = as.character(Movement_category),
         Movement_category = case_when(
           Movement_category == "Linear Distance - Speed" ~ "Linear Distance",
               Movement_category == "Tortuosity - Opposite Direction" ~ "Tortuosity",
           TRUE ~ Movement_category
         ),
         Movement_category = as.factor(Movement_category))

table(data$Movement_category)    

## Cleaning Fire Category

# looking at unique category counts to identify duplicates 

table(data$Impact_fire_type)

data$Impact_fire_type <- factor(data$Impact_fire_type)

data$class_ncbi = factor(data$class_ncbi)

table(data$class_ncbi)

data <- data %>%
  mutate(class_ncbi = if_else(class_ncbi == "Lepidosauria",
                              "Reptilia",
                              class_ncbi))

# Changing NAs from the torts to reptiles

data$class_ncbi[is.na(data$class_ncbi)] <- "Reptilia"

## Writing CSV

write_csv(data, "Outputs/Dat_for_MA_analysisRev5.csv")

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

print(es) 

# obs level random effect

es$Obs_ID <- factor(1:nrow(es))

# Study level random effect

es$StudyID <- NULL

es$Study_ID <- factor(as.numeric(as.factor(es$Firstauthor_Year)))

# Correcting direction of Tortuouscity Harris2020 

es %>%
  dplyr::filter(Firstauthor_Year %in% c("Harris2020"))
  
es <- es %>%
  mutate(yi = ifelse(Obs_ID %in% c(91), yi * -1, yi))


################################## MA Models ##################################

# Variance covariance matrix

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5,
             data = es) 

#Overall model

rerun <- FALSE

if(rerun){
  mod.overall <- rma.mv(yi = yi, V = vcv,
                        random = list(~1 | Study_ID,
                                      ~1 | Species_tree, # phylo effect 
                                      ~1 | Species2, # non-phylo effect 
                                      ~1 | Obs_ID), 
                        data =  es,
                        control = list(optimizer="BFGS"),
                        test = "t",
                        sparse = TRUE,
                        R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall, here("Outputs/No_space_use", "mod.overall.Rds"))
  summary(mod.overall)
}else{
  mod.overall <- readRDS(here("Outputs/No_space_use", "mod.overall.Rds"))
}

I2 = round(i2_ml(mod.overall), 2)

overall <- orchard_plot(mod.overall, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                        angle = 0)
overall

ggsave(here("Outputs/No_space_use", "Overall_nospace.pdf"), width = 6, height = 4)

# Heterogeneity of the effect sizes 

round(i2_ml(mod.overall), 2)

###############################################################################
###################### Moderator: Movement Category ###########################
###############################################################################

if(rerun){
  mod.overall_move <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Study_ID,
                                           ~1 | Species_tree, # phylo effect 
                                           ~1 | Species2, # non-phylo effect 
                                           ~1 | Obs_ID), 
                             data =  es,
                             mods = ~ Movement_category -1,
                             control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall_move, here("Outputs/No_space_use", "mod.overall_move.Rds"))
}else{
  mod.overall_move <- readRDS(here("Outputs/No_space_use", "mod.overall_move.Rds"))
}

summary(mod.overall_move)

# Plotting 

overall_move <- orchard_plot(mod.overall_move, xlab = "Change in movement (Hedge's g)", group = "Study_ID", mod = "Movement_category",
                             angle = 0)
overall_move

# Saving Plot

ggsave(here("Outputs/No_space_use", "Overall_move_nospace.pdf"), width = 6, height = 4)

###############################################################################
###################### Moderator: Species #####################################
###############################################################################

if(rerun){
  mod.species <- rma.mv(yi = yi, V = vcv,
                        random = list(~1 | Study_ID,
                                      ~1 | Species_tree, # phylo effect 
                                      ~1 | Species2, # non-phylo effect 
                                      ~1 | Obs_ID), 
                        data =  es,
                        mods = ~ Species2 -1,
                        control = list(optimizer="BFGS"),
                        test = "t",
                        sparse = TRUE,
                        R = list(Species_tree = cor1)
  )
  saveRDS(mod.species, here("Outputs/No_space_use", "mod_species.Rds"))
}else{
  mod.species <- readRDS(here("Outputs/No_space_use", "mod_species.Rds"))
}

summary(mod.species)

# Plotting 

# Creating Palette 
n_species <- 40  
n_species
random_palette <- sample(viridis(n_species, option = "D"))

#Plot

Species <- suppressMessages(
  orchard_plot(
    mod.species,
    xlab = "Change in movement (Hedge's g)",
    group = "Study_ID",
    mod = "Species2",
    angle = 0
  ) +
    scale_fill_manual(values = random_palette) +
    scale_color_manual(values = random_palette) +
    theme(axis.text.y = element_text(hjust = 0)) +
    theme(legend.position = "bottom", legend.justification = "right")
)

Species

ggsave(here("Outputs", "No_space_use", "Species.pdf"), width = 14, height = 12)


###############################################################################
############################# Moderator: Fire Type ############################
###############################################################################

# Filtering out not provided, unclear or both fire types

fire.es <- es %>%
  filter(Impact_fire_type %in% c("Planned", "Wildfire"))

fire.es$Impact_fire_type <- factor(fire.es$Impact_fire_type)

# Variance covariance matrix

fire.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5,
                  data = fire.es) 

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
                     R = list(Species_tree = cor1)
  )
  saveRDS(mod.fire, here("Outputs/No_space_use", "mod_fire.Rds"))
}else{
  mod.fire <- readRDS(here("Outputs/No_space_use", "mod_fire.Rds"))
}

summary(mod.fire)

# Plotting 

plot.fire <- orchard_plot(mod.fire, xlab = "Change in movement (Hedge's g)", 
                          group = "Study_ID", mod = "Impact_fire_type", angle = 0)
plot.fire

# Saving Plot

ggsave(here("Outputs", "No_space_use", "Fire_type.pdf"), width = 6, height = 4)


###############################################################################
######################## Moderator: Class #####################################
###############################################################################

table(es$class_ncbi)

# Filtering out Amphibia due to small sample size 

class.es <- es %>%
  filter(class_ncbi %in% c("Aves", "Mammalia", "Reptilia"))

class.es$class_ncbi <- factor(class.es$class_ncbi)

# VCV

class.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                   data = class.es) 

# Model
if(rerun){
  mod.class <- rma.mv(yi = yi, V = class.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  class.es,
                      mods = ~ class_ncbi -1, 
                      control = list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1)
  )
  saveRDS(mod.class, here("Outputs/No_space_use", "mod_class.Rds"))
}else{
  mod.class <- readRDS(here("Outputs/No_space_use", "mod_class.Rds"))
}

  summary(mod.class)

plot.class <- orchard_plot(mod.class, xlab = "Change in movement (Hedge's g)", 
                           group = "Study_ID", mod = "class_ncbi", angle = 0)
plot.class

ggsave(here("Outputs", "No_space_use", "class.pdf"), width = 6, height = 4)

###############################################################################
######################### Biome as Moderator ##################################
###############################################################################

table(es$biome_name)

# Filtering out biomes with <10 es

biome.es <- es %>%
  group_by(biome_name) %>%
  filter(n() >= 10) %>%
  ungroup()

table(biome.es$biome_name)

biome.es$biome_name <- factor(biome.es$biome_name)

# VCV

biome.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                   data = biome.es) 

# Model

if(rerun){
  mod.biome <- rma.mv(yi = yi, V = biome.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  biome.es,
                      mods = ~ biome_name -1, #' [We remove the intercept here bc we want estimates for all values and not a comparison to the first one]
                      control=list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1) #' [EW - this is the same as above]
  )
  saveRDS(mod.biome, here("Outputs/No_space_use", "mod_biome.Rds"))
}else{
  mod.biome <- readRDS(here("Outputs/No_space_use", "mod_biome.Rds"))
}

summary(mod.biome)

# Plotting

plot.biome <- orchard_plot(mod.biome, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                           mod = "biome_name", angle = 0)

plot.biome

#Fixing labels

wrapped_labels <- function(x) str_wrap(x, width = 25)

plot.biome + 
  scale_x_discrete(labels = wrapped_labels) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 8)) +
  theme(axis.title.x = element_text(size = 10))

# Saving Plot

ggsave(here("Outputs", "No_space_use", "Biome.pdf"), width = 8, height = 6)

###############################################################################
######################### Veg_clim as Moderator ###############################
###############################################################################

table(es$veg_climate)

# Filtering out closed_Arid and open_Temperate due to small sample size (<10)

veg_clim.es <- es %>%
  filter(veg_climate %in% c("closed_Cold", 
                            "closed_Temperate",
                            "open_Arid", 
                            "open_Tropical" ))


veg_clim.es$veg_climate <- factor(veg_clim.es$veg_climate)

# VCV

veg_clim.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                      data = veg_clim.es) 

# Model

if(rerun){
  mod.veg_clim <- rma.mv(yi = yi, V = veg_clim.vcv,
                         random = list(~1 | Study_ID,
                                       ~1 | Species_tree, # phylo effect 
                                       ~1 | Species2, # non-phylo effect 
                                       ~1 | Obs_ID), 
                         data =  veg_clim.es,
                         mods = ~ veg_climate -1, #' 
                         control=list(optimizer="BFGS"),
                         test = "t",
                         sparse = TRUE,
                         R = list(Species_tree = cor1) 
  )
  saveRDS(mod.veg_clim, here("Outputs/No_space_use", "mod_veg_clim.Rds"))
}else{
  mod.veg_clim <- readRDS(here("Outputs/No_space_use", "mod_veg_clim.Rds"))
}

summary(mod.veg_clim)

# Plotting

plot.veg_clim <- orchard_plot(mod.veg_clim, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                              mod = "veg_climate", angle = 0)

plot.veg_clim

# Saving Plot

ggsave(here("Outputs", "No_space_use", "Veg_clim.pdf"), width = 8, height = 6)

###############################################################################
######################### Climate as Moderator ###############################
###############################################################################

table(es$Climate)

es$Climate <- factor(es$Climate)

# Model

if(rerun){
  mod.clim <- rma.mv(yi = yi, V = vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  es,
                     mods = ~ Climate -1, #' [We remove the intercept here bc we want estimates for all values and not a comparison to the first one]
                     control=list(optimizer="BFGS"),
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1) 
  )
  saveRDS(mod.clim, here("Outputs/No_space_use", "mod_clim.Rds"))
}else{
  mod.clim <- readRDS(here("Outputs/No_space_use", "mod_clim.Rds")) #' error here - saving wrong model - fixed
}
summary(mod.clim)

# Plotting

plot.clim <- orchard_plot(mod.clim, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                          mod = "Climate", angle = 0)

plot.clim

# Saving Plot

ggsave(here("Outputs", "No_space_use", "clim.pdf"), width = 8, height = 6)

###############################################################################
######################## Body Mass model : Overall  ###########################
###############################################################################

# Tranforming data

es$logmass <- log(es$BodyMass.g)

if(rerun){
  mod.bm <- rma.mv(yi = yi, V = vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  es,
                   mods = ~ logmass,
                   control=list(optimizer="BFGS"),
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm, here("Outputs/No_space_use", "mod1_bm.Rds"))
}else{
  mod.bm <- readRDS(here("Outputs/No_space_use", "mod1_bm.Rds"))
}

summary(mod.bm)

#Extracting and plotting body mass models 

round(r2_ml(mod.bm)*100, 2) #This is r2 of the line

## es[is the data frame]$logmass [is the body mass column] 

newdat <- data.frame(logmass = seq(min(es$logmass, na.rm = TRUE),
                                   max(es$logmass, na.rm = TRUE),
                                   length.out = 100))

X <- model.matrix(~ logmass, data = newdat)[, -1, drop = FALSE]

preds <- predict(mod.bm, newmods = X)

plot_data <- cbind(newdat, preds)

setDT(es)

es[, wi := 1/sqrt(vi)]
es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
es

prey_bm_text <- tidy(mod.bm)
prey_bm_text

# Plotting 

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

ggsave(here("Outputs/No_space_use", "body mass.pdf"), height = 4, width = 8)


###############################################################################
######################## Body Mass model : Mammals  ###########################
###############################################################################

bm.mam.es <- es %>%
  filter(class_ncbi %in% c("Mammalia"))

#VCV

mam.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                 data = bm.mam.es) 

# Model 

if(rerun){
  mod.bm.mam <- rma.mv(yi = yi, V = mam.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.mam.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.mam, here("Outputs/No_space_use", "mod_bm_mam.Rds"))
}else{
  mod.bm.mam <- readRDS(here("Outputs/No_space_use", "mod_bm_mam.Rds"))
}

summary(mod.bm.mam)

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.mam)*100, 2) #This is r2 of the line

## es[is the data frame]$logmass [is the body mass column] 

newdat.mam <- data.frame(logmass = seq(min(bm.mam.es$logmass, na.rm = TRUE),
                                       max(bm.mam.es$logmass, na.rm = TRUE),
                                       length.out = 100)) 

X <- model.matrix(~ logmass, data = newdat.mam)[, -1, drop = FALSE]

preds.mam <- predict(mod.bm.mam, newmods = X)

plot_data.mam <- cbind(newdat.mam, preds.mam)

setDT(bm.mam.es)

bm.mam.es[, wi := 1/sqrt(vi)]
bm.mam.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.mam.es

prey_bm.mam_text <- tidy(mod.bm.mam)
prey_bm.mam_text

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.mam.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.mam.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.mam,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.mam, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.mam, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Mammal Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.mam.mean

ggsave(here("Outputs/No_space_use", "body mass mammal.pdf"), height = 4, width = 8)

###############################################################################
######################## Body Mass model : Birds  #############################
###############################################################################

bm.ave.es <- es %>%
  filter(class_ncbi %in% c("Aves"))

# VCV
aves.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5,
                  data = bm.ave.es) 
# Model

if(rerun){
  mod.bm.ave <- rma.mv(yi = yi, V = aves.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.ave.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.ave, here("Outputs/No_space_use", "mod_bm_ave.Rds"))
}else{
  mod.bm.ave <- readRDS(here("Outputs/No_space_use", "mod_bm_ave.Rds"))
}

summary(mod.bm.ave)

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.ave)*100, 2) #This is r2 of the line

## es[is the data frame]$logmass [is the body mass column] 

newdat.ave <- data.frame(logmass = seq(min(bm.ave.es$logmass, na.rm = TRUE),
                                       max(bm.ave.es$logmass, na.rm = TRUE),
                                       length.out = 100)) 

X <- model.matrix(~ logmass, data = newdat.ave)[, -1, drop = FALSE]

preds.ave <- predict(mod.bm.ave, newmods = X)

plot_data.ave <- cbind(newdat.ave, preds.ave)

setDT(bm.ave.es)

bm.ave.es[, wi := 1/sqrt(vi)]
bm.ave.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.ave.es

prey_bm.ave_text <- tidy(mod.bm.ave)
prey_bm.ave_text

# Plotting

p.bm.ave.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.ave.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.ave,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.ave, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.ave, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Aves Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.ave.mean

ggsave(here("Outputs/No_space_use", "body mass aves.pdf"), height = 4, width = 8)

###############################################################################
######################## Body Mass model : Reptiles ###########################
###############################################################################

bm.rep.es <- es %>%
  filter(class_ncbi %in% c("Reptilia"))

# VCV

rep.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                 data = bm.rep.es) 
# Model

if(rerun){
  mod.bm.rep <- rma.mv(yi = yi, V = rep.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.rep.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.rep, here("Outputs/No_space_use", "mod_bm_rep.Rds"))
}else{
  mod.bm.rep <- readRDS(here("Outputs/No_space_use", "mod_bm_rep.Rds"))
}

summary(mod.bm.rep)

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.rep)*100, 2) #This is r2 of the line

## es[is the data frame]$logmass [is the body mass column] 

newdat.rep <- data.frame(logmass = seq(min(bm.rep.es$logmass, na.rm = TRUE),
                                       max(bm.rep.es$logmass, na.rm = TRUE),
                                       length.out = 100)) 

X <- model.matrix(~ logmass, data = newdat.rep)[, -1, drop = FALSE]

preds.rep <- predict(mod.bm.rep, newmods = X)

plot_data.rep <- cbind(newdat.rep, preds.rep)

setDT(bm.rep.es)

bm.rep.es[, wi := 1/sqrt(vi)]
bm.rep.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.rep.es

prey_bm.rep_text <- tidy(mod.bm.rep)
prey_bm.rep_text 

# plotting

p.bm.rep.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.rep.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.rep,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.rep, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.rep, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Reptilia Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.rep.mean

ggsave(here("Outputs/No_space_use", "body mass rept.pdf"), height = 4, width = 8)

###############################################################################
############################ Burn Treatment Age Model #########################
###############################################################################

# Coercing variable to integer,looking for NA's
es <- es %>%
  mutate(Maximum_treatment_burn_age_days = if_else(
    grepl("^[-]?[0-9]+$", Maximum_treatment_burn_age_days),
    as.integer(Maximum_treatment_burn_age_days),
    NA_integer_
  ))

# Looking for NAs

summary(es$Maximum_treatment_burn_age_days)

# Removing NAs

age.es <- es[!is.na(es$Maximum_treatment_burn_age_days), ]

age.es <- es %>%
  mutate(Maximum_treatment_burn_age_days = if_else(
    Maximum_treatment_burn_age_days == 60989, 
    609, 
    Maximum_treatment_burn_age_days
  ))

summary(age.es$Maximum_treatment_burn_age_days)

hist(age.es$Maximum_treatment_burn_age_days)

age.es$log.age <- log(age.es$Maximum_treatment_burn_age_days)

hist(age.es$log.age)

# Model

if(rerun){
  mod.age <- rma.mv(yi = yi, V = vcv,
                    random = list(~1 | Study_ID,
                                  ~1 | Species_tree, # phylo effect 
                                  ~1 | Species2, # non-phylo effect 
                                  ~1 | Obs_ID), 
                    data =  age.es,
                    mods = ~ log.age, #uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                    test = "t",
                    sparse = TRUE,
                    R = list(Species_tree = cor1)
  )
  saveRDS(mod.age, here("Outputs/No_space_use", "mod1_age.Rds"))
}else{
  mod.age <- readRDS(here("Outputs/No_space_use", "mod1_age.Rds"))
}

summary(mod.age)

#Extracting and plotting

round(r2_ml(mod.age)*100, 2) #This is r2 of the line

## es[is the data frame]$logmass [is the body mass column] 

newdat.age <- data.frame(log.age = seq(min(age.es$log.age, na.rm = TRUE),
                                       max(age.es$log.age, na.rm = TRUE),
                                       length.out = 100))

X.age <- model.matrix(~ log.age, data = newdat.age)[, -1, drop = FALSE]

preds.age <- predict(mod.age, newmods = X.age)

plot_data.age <- cbind(newdat.age, preds.age)

setDT(age.es)

age.es[, wi := 1/sqrt(vi)]
age.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
age.es

prey_age_text <- tidy(mod.age)
prey_age_text

p.age.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = age.es, 
              aes(x = log.age, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.age,
              aes(x = log.age,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.age, 
              aes(x = log.age,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.age, 
            aes(x = log.age,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log Maximum Burn Age (days)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.age.mean #' [EW - i get Removed 2 rows containing missing values or values outside the scale range (`geom_point()`). need to reset the rang of values i think although it's bc you don't have coords set. Not a big deal anyway.]

ggsave(here("Outputs/No_space_use", "Burn Age.pdf"), height = 4, width = 8)

###############################################################################
########################## Model: Fire Size ###################################
###############################################################################

# Coercing variable to integer, looking for NA's

size.es <- es %>%
  mutate(Impact_fire_size_ha = as.numeric(Impact_fire_size_ha) ) %>%
  mutate(Impact_fire_size_ha = as.integer(round(Impact_fire_size_ha)))

# Looking for NAs
summary(size.es$Impact_fire_size_ha)

# Removing NAs
size.es <- size.es[!is.na(size.es$Impact_fire_size_ha), ]

# Transforming data

hist(size.es$Impact_fire_size_ha)

size.es$log.size <- log(size.es$Impact_fire_size_ha)

hist(size.es$log.size)

#VCV 

size.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, 
                  data = size.es) 
#Model

if(rerun){
  mod.size <- rma.mv(yi = yi, V = size.vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  size.es,
                     mods = ~ log.size, 
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1)
  )
  saveRDS(mod.size, here("Outputs/No_space_use", "mod1_size.Rds"))
}else{
  mod.size <- readRDS(here("Outputs/No_space_use", "mod1_size.Rds"))
}

summary(mod.size)

round(r2_ml(mod.size)*100, 2) #This is r2 of the line

# Plotting

## es[is the data frame]$logmass [is the body mass column] 

newdat.size <- data.frame(log.size = seq(min(size.es$log.size, na.rm = TRUE),
                                         max(size.es$log.size, na.rm = TRUE),
                                         length.out = 100)) 

X.size <- model.matrix(~ log.size, data = newdat.size)[, -1, drop = FALSE]

preds.size <- predict(mod.size, newmods = X.size)

plot_data.size <- cbind(newdat.size, preds.size)

setDT(size.es)

size.es[, wi := 1/sqrt(vi)]
size.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
size.es

prey_age_text <- tidy(mod.size)
prey_age_text

p.size.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = size.es, 
              aes(x = log.size, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.size,
              aes(x = log.size,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.size, 
              aes(x = log.size,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.size, 
            aes(x = log.size,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log Burn Area (ha)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.size.mean

ggsave(here("Outputs/No_space_use", "Burn Size.pdf"), height = 4, width = 8)

################################################################################
################### Pausas2017 Fire Activity ###################################
################################################################################

summary(es$fireactivi)

hist(es$fireactivi) #trimodal, may need to run as categorical rather than continuous if we use?

fireact.es <- es

# Model

if(rerun){
  mod.fireact <- rma.mv(yi = yi, V = vcv,
                        random = list(~1 | Study_ID,
                                      ~1 | Species_tree, # phylo effect 
                                      ~1 | Species2, # non-phylo effect 
                                      ~1 | Obs_ID), 
                        data =  es,
                        mods = ~ fireactivi,
                        control=list(optimizer="BFGS"),
                        test = "t",
                        sparse = TRUE,
                        R = list(Species_tree = cor1)
  )
  saveRDS(mod.fireact, here("Outputs/No_space_use", "mod_fireact.Rds"))
  summary(mod.fireact)
}else{
  mod.fireact <- readRDS(here("Outputs/No_space_use", "mod_fireact.Rds"))
}

mod.fireact

round(r2_ml(mod.fireact)*100, 2) #This is r2 of the line

# Plotting

newdat.fireact <- data.frame(fireactivi = seq(min(es$fireactivi, na.rm = TRUE), ## [the data frame]$[the variable] 
                                         max(es$fireactivi, na.rm = TRUE),
                                         length.out = 100)) 

X.fireactivity <- model.matrix(~ fireactivi, data = newdat.fireact)[, -1, drop = FALSE]

preds.fireact <- predict(mod.fireact, newmods = X.fireactivity)

plot_data.fireact <- cbind(newdat.fireact, preds.fireact)

setDT(es)

es[, wi := 1/sqrt(vi)]
es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
es

prey_age_text <- tidy(mod.fireact)
prey_age_text

p.fireact.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = es, 
              aes(x = fireactivi, size = pt_size, # do I need to change? 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.fireact,
              aes(x = fireactivi,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.fireact, 
              aes(x = fireactivi,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.fireact, 
            aes(x = fireactivi,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log Fire Activity Index")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.fireact.mean

ggsave(here("Outputs/No_space_use", "Fire activity.pdf"), height = 4, width = 8)

################################################################################
########################### Publication bias analyses #########################
################################################################################

es$effectN <- (es$Unburned_n * es$Burned_n) / (es$Unburned_n + es$Burned_n)
es$sqeffectN <- sqrt(es$effectN)

if(rerun){
  mod_diff_egg <- rma.mv(yi = yi, V = vcv,
                         random = list(~1 | Study_ID,
                                       ~1 | Species_tree, # phylo effect 
                                       ~1 | Species2, # non-phylo effect 
                                       ~1 | Obs_ID), 
                         data =  es,
                         mods = ~ sqeffectN,
                         control=list(optimizer="BFGS"),
                         test = "t",
                         sparse = TRUE,
                         R = list(Species_tree = cor1)
  )
  summary(mod_diff_egg)
  saveRDS(mod_diff_egg, here("Outputs/No_space_use", "eggers.Rds"))
}else{
  mod_diff_egg <- readRDS(here("Outputs/No_space_use", "eggers.Rds"))
}

mod_diff_egg #No significant intercept

##Funnel

funnel(mod_diff_egg, 
       yaxis="seinv",
       xlab = "Standardized residuals",
       ylab = "Precision (inverse of SE)",
       
)

ggsave("Outputs/No_Space_use/funnel_main.pdf", width = 5, height = 5)

### Decline effect for difference

str(es$Publication_Year)

mod_mad <- rma.mv(yi = yi, V = vcv,
                  random = list(~1 | Study_ID,
                                ~1 | Species_tree, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data = es,
                  mods = ~ Publication_Year,
                  control=list(optimizer="BFGS"),
                  test = "t",
                  sparse = TRUE,
                  R = list(Species_tree = cor1)
)

summary(mod_mad)
saveRDS(mod_mad, here("Outputs/No_space_use", "decline.Rds"))

# Plotting

decline <- bubble_plot(mod_mad,
                       mod = "Publication_Year",
                       group = "Study_ID",
                       xlab = "Publication year",
                       ylab = "Hedge's g",
                       g = TRUE,
                       weights = 1/sqrt(es$vi))

decline ## if this is giving an error, just restart R and it should work

ggsave("Outputs/No_space_use/decline.pdf", width = 8, height = 4)


### >>>> leave-one-out analysis 

##### if rerunning load the model below ####

es$Firstauthor_Year <- as.factor(es$Firstauthor_Year)

es$leave_out <- as.factor(es$Firstauthor_Year)

LeaveOneOut_effectsize <- list()
for (i in 1:length(levels(es$leave_out))) {
  temp_dat <- es %>%
    dplyr::filter(leave_out != levels(es$leave_out)[i])
  
  VCV_leaveout <- vcalc(vi = temp_dat$vi, cluster = temp_dat$Study_ID, obs = temp_dat$Obs_ID, rho = 0.5)
  
  LeaveOneOut_effectsize[[i]] <-  rma.mv(yi = yi, V = VCV_leaveout,
                                         random = list(~1 | Study_ID,
                                                       ~1 | Species_tree, # phylo effect 
                                                       ~1 | Species2, # non-phylo effect 
                                                       ~1 | Obs_ID), 
                                         control=list(optimizer="BFGS"),
                                         test = "t",
                                         method = "REML", 
                                         sparse = TRUE,
                                         R = list(Species_tree = cor1),
                                         data = temp_dat[temp_dat$leave_out != levels(temp_dat$leave_out)[i], ])
}


# writing function for extracting est, ci.lb, and ci.ub from all models
est.func <- function(model) {
  df <- data.frame(est = model$b, lower = model$ci.lb, upper = model$ci.ub)
  return(df)
}

# using dplyr to form data frame
MA_oneout <- lapply(LeaveOneOut_effectsize,function(x) est.func(x)) %>%
  bind_rows %>%
  mutate(left_out = levels(es$leave_out))

# telling ggplot to stop reordering factors
MA_oneout$left_out <- as.factor(MA_oneout$left_out)
MA_oneout$left_out <- factor(MA_oneout$left_out, levels = MA_oneout$left_out)

# saving the runs
saveRDS(MA_oneout, here("Outputs/No_space_use", "MA_oneout.RDS"))

#Load model here if rerunning 

MA_oneout <- readRDS(here("Outputs/No_space_use", "MA_oneout.RDS"))
summary(MA_oneout)

# plotting
leaveoneout <- ggplot(MA_oneout) + geom_hline(yintercept = 0, lty = 2, lwd = 1) +
  geom_hline(yintercept = mod.overall$ci.lb, lty = 3, lwd = 0.75, colour = "black") +
  geom_hline(yintercept = mod.overall$b, lty = 1, lwd = 0.75, colour = "black") + 
  geom_hline(yintercept = mod.overall$ci.ub,
             lty = 3, lwd = 0.75, colour = "black") + 
  geom_pointrange(aes(x = left_out, y = est,
                      ymin = lower, ymax = upper)) + 
  xlab("Study left out") + 
  ylab("Difference in logit overlap (logit)") +
  coord_flip() + 
  theme(panel.grid.minor = element_blank()) + theme_bw() + theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) + theme(axis.text.y = element_text(size = 6))

leaveoneout

ggsave("Outputs/No_space_use/leaveoneout.pdf", width = 8, height = 4)

#' [END]

#for making tables for the supplement

mod1 <- coef(summary(mod.overall))#paste model names into the brackets
mod2 <- coef(summary(mod.overall_move))
mod3 <- coef(summary(mod.class))
 
mod4 <- coef(summary(mod.fire))
mod5 <- coef(summary(mod.size))
mod6 <- coef(summary(mod.age))

mod7 <- coef(summary(mod.biome))

mod8 = coef(summary(mod.bm))
mod9 <- coef(summary(mod.bm.mam))
mod10 <- coef(summary(mod.bm.ave))
mod11 <- coef(summary(mod.bm.rep))

mod12 = coef(summary(mod.veg_clim))
mod13 = coef(summary(mod.clim))

mod14 = coef(summary(mod.fireact))

mod15 = coef(summary(mod_diff_egg))
mod16 = coef(summary(mod_mad))


table1 <- rbind(mod1,mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16) #bind the models together 

table1$term <- rownames(table1)

rownames(table1) <- NULL

write.csv(table1, "Outputs/No_space_use/Sup_tableRev5.csv") # Write it out to CSV to then paste into word doc


