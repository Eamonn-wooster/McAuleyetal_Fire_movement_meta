##### Superseded ############################################################## 

## Joining sample sizes back into condensed dataset ##



## CM

## General Packages ##

library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)


## Setting directory ##

here()

## Reading in data with sample sizes
dat= read.csv(here("Data/Extraction_data_uncondensed.csv"))

## Reading in condensed data

dat.condensed = read.csv(here("Data/Studies_Converted_Condensed-Rev2.csv"))

## Joining sample sizes by Study ID

dat.condensed.n <- dat.condensed %>%
  left_join(dat %>% select(Study_ID, Sample_Size_Before, Sample_Size_After, 
                           Sample_Size_Control, Sample_Size_Impact, Sample_Size_Control_Before,
                           Sample_Size_Control_After, Sample_Size_Impact_Before, Sample_Size_Impact_After),
            by = "Study_ID")
## Merging CB and IA into burned and unburned sample sizes

dat.condensed.n <- dat.condensed.n %>%
  mutate(Burned_Sample_Size = coalesce(Sample_Size_After, Sample_Size_Impact)) %>%
  mutate(Unburned_Sample_Size = coalesce(Sample_Size_Before, Sample_Size_Control))

## Summing sample sizes for BACI studies into Burned and unburned 

dat.condensed.n <- dat.condensed.n %>%
  mutate(BACI_Sample_Size_Unburned = Sample_Size_Control_Before + Sample_Size_Control_After) %>%
  mutate(BACI_Sample_Size_Burned = Sample_Size_Impact_Before + Sample_Size_Impact_After)

## Merging BACI, CB and AI sample sizes into Burned and Unburned sample sizes

dat.condensed.n = dat.condensed.n %>%
  mutate(All_Burned_Sample_Size = coalesce(Burned_Sample_Size, BACI_Sample_Size_Burned)) %>%
  mutate(All_Unburned_Sample_Size = coalesce(Unburned_Sample_Size, BACI_Sample_Size_Unburned))

## Removing unnecessary columns

dat.condensed.n <- dat.condensed.n %>%
  select(-Sample_Size_Control_Before, -Sample_Size_Control_After, -Sample_Size_Impact_Before,
         -Sample_Size_Impact_After, -Sample_Size_Before, -Sample_Size_After, -Sample_Size_Control, -Sample_Size_Impact)

dat.condensed.n <- dat.condensed.n %>%
  select(-BACI_Sample_Size_Unburned, -BACI_Sample_Size_Burned)

dat.condensed.n <- dat.condensed.n %>%
  select(-Burned_Sample_Size, -Unburned_Sample_Size)

## Renaming sample size columns

dat.condensed.n <- dat.condensed.n %>%
  rename(Burned_n = All_Burned_Sample_Size) %>%
  rename(Unburned_n = All_Unburned_Sample_Size)

## Checking for NAs in sample size

dat.condensed.n %>%
  filter(is.na(Burned_n) | is.na(Unburned_n))

## Writing csv file

write_csv(dat.condensed.n, "Outputs/Dat_for_analysis-Rev2.csv")

## Writing current data_for_analysis.csv file - 26/04/2025

write.csv(dat.condensed.n, "Outputs/Dat_for_analysis_CM.csv")
