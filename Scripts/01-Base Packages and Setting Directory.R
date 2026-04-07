## Packages ##

library(tidyverse)
library(here)
library(dplyr)

## Setting directory ##

here()

## Reading in data with sample sizes

dat= read.csv(here("Data/Data_Extracted_May-2025-MASTER.csv"))

## Reading in condensed data

dat.condensed = read.csv(here("Data/Studies_Converted_Condensed-Rev3.csv"))

## Joining sample sizes by Study ID

dat.condensed.n <- dat.condensed %>%
  left_join(dat %>% select(Study_ID, Sample_Size_Before, Sample_Size_After, 
            Sample_Size_Control, Sample_Size_Impact, Sample_Size_Control_Before, Sample_Size_Impact_Before),
            by = "Study_ID")

# Merging C-B and I-A into burned and unburned sample sizes (Following Morris 2008, 
# and using before sample sizes and pooled SD to calculate effect sizes with hedges g)

dat.condensed.n <- dat.condensed.n %>%
  mutate(Burned_Sample_Size = coalesce(Sample_Size_After, Sample_Size_Impact, Sample_Size_Impact_Before)) %>%
  mutate(Unburned_Sample_Size = coalesce(Sample_Size_Before, Sample_Size_Control, Sample_Size_Impact_Before))

############ Not doing this anymore ###########################################

## Summing sample sizes for BACI studies into Burned and unburned 

#dat.condensed.n <- dat.condensed.n %>%
  #mutate(BACI_Sample_Size_Unburned = Sample_Size_Control_Before + Sample_Size_Control_After) %>%
  #mutate(BACI_Sample_Size_Burned = Sample_Size_Impact_Before + Sample_Size_Impact_After)

## Merging BACI, CB and AI sample sizes into Burned and Unburned sample sizes

## dat.condensed.n = dat.condensed.n %>%
  #mutate(All_Burned_Sample_Size = coalesce(Burned_Sample_Size, BACI_Sample_Size_Burned)) %>%
  #mutate(All_Unburned_Sample_Size = coalesce(Unburned_Sample_Size, BACI_Sample_Size_Unburned))

############# End of no longer doing this #######################################################

## Removing unnecessary columns

dat.condensed.n <- dat.condensed.n %>%
  select(-Sample_Size_Control_Before, -Sample_Size_Impact_Before,
         , -Sample_Size_Before, -Sample_Size_After, -Sample_Size_Control, -Sample_Size_Impact)

## Renaming sample size columns

dat.condensed.n <- dat.condensed.n %>%
  rename(Burned_n = Burned_Sample_Size) %>%
  rename(Unburned_n = Unburned_Sample_Size)

## Checking NAs and sample sizes

dat.condensed.n %>%
  filter(is.na(Burned_n) | is.na(Unburned_n))

any(dat.condensed.n$Burned_n < 5, na.rm = TRUE)

any(dat.condensed.n$Unburned_n < 5, na.rm = TRUE)

dat.condensed.n %>%
  filter(Burned_n < 5)

dat.condensed.n %>%
  filter(Unburned_n <5)

## Blakey2022 included due to resampling method
## Shaw1990 included as 46 bison tagged, but all grouped together for analysis, so n is number of months observed 
#(ie proportion of observations inside and outside not provided for each individual)

## Removing Kortner2007 and Michaud2017 rows with <5

dat.condensed.n <- dat.condensed.n %>%
  filter(!Study_ID %in% c("165", "5"))

## Writing csv file

write_csv(dat.condensed.n, "Outputs/Dat_for_analysis-ss-Rev3.csv")
