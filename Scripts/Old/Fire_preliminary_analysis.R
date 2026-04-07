# pilot analysis - updating with new col names for revision
rm(list = ls())

# install.packages("pacman")

# devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
# pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans,
#               ape, phytools, flextable)

### Loading packages ----
pacman::p_load(devtools, 
               tidyverse, 
               metafor, 
               patchwork, 
               R.rsp, 
               emmeans,
               orchaRd,
               metafor,
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

#getwd()

# using here to make this code work for both Mac and Windowns
#This should work for you now theres a proj? 

here()

dat <- read.csv(here("Outputs/Dat_for_analysis-ss-aves-rep-mass-Rev3.csv")) 

# getting tree

tree1 <- read.tree(here("Outputs/Tree_fire-Rev3.tre")) 


#These tree are fine - making new trees for each evidence catergories has no influence on the model
# getting branch length and correlation matrix
tree1b <- compute.brlen(tree1)
cor1 <-  vcv(tree1b, corr=T)


# creating predator and pray names which match with the trees

# # checking the match


# # checking the match - matches : ) 
setdiff(dat$Species_tree, tree1$tip.label)


# # creating non-phylo columns
dat$Species2 <- dat$Species_tree

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

#Taking the upper value from the mean gives us the SD values
SE$sd_Unburned <- SE$Upper_Unburned - SE$Mean_Unburned

SE$sd_Fire <- SE$Upper_Fire - SE$Mean_Fire

#now convert SE to SD

SE$sd_Unburned <- se_to_sd(SE$sd_Unburned, SE$Unburned_n)

SE$sd_Fire <- se_to_sd(SE$sd_Fire, SE$Burned_n)

#Put them back together!! 

data <- rbind(SD, CI, SE)

### CM Hedges g effect size, with help from chatgpt ###

##CM there are negative SDs to be sorted out... 
##CM removing these
data_clean <- subset(data, sd_Unburned > 0 & sd_Fire > 0)

# Calculate Hedges' g
es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = data_clean,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results
print(es) ## CM yi is the SMD effect size, and vi is the sampling variance. The variance here seems better??

####### End chatgpt help #####

# #ok! effect size
# 
data$RR <- log(data$Mean_Fire/data$Mean_Unburned)
# 
# This is following Doherty - maybe check i've pulled this formula across properly!! 
data$vi <- ((data$sd_Fire^2)/data$Unburned_n + data$Mean_Fire^2) + 
  ((data$sd_Unburned^2)/data$Unburned_n + data$Mean_Unburned^2)

# one needs to transform variance too using the delta method - do this? Yes, do this and update Eamonn

data$vi_d <- data$vi  / (data$RR * (1 - data$RR))^2



##### Models--------


# >>> Overall model ----------------------------------------------

vcv2 <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
              data = data) 

#for the model loops 
rerun <- F
if(rerun){
  mod1 <- rma.mv(yi = overlap_diff, V = vcv2,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect of Predactor
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), # this is needed
                      data =  dat,
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1)
  )
  saveRDS(mod1_diff, here("R", "Outputs", "mod1_overall_difference.Rds"))
  summary(mod1_diff)
}else{
  mod1_diff <- readRDS(here("R", "Outputs", "mod1_overall_difference.Rds"))
}

dif <- orchard_plot(mod1_diff, xlab = "dff in logit(overlap %)", group = "Study_ID")

dif

ggsave(here("Figures", "overlap_overall.pdf"), width = 6, height = 4)

