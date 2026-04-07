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
               data.table,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp #,#added
               # miWQS #added
) 

# Packages #
remotes::install_github("RS-eco/geodat", force = TRUE)
library(geodat)
library(sf)     
library(dplyr) 
library(tidyr)    
library(stringr)  
library(here)


# Setting WD

here()

## Importing extraction data

LR.dat = read.csv(here("Data/Literature Review/LR_Data-June-2025.csv"))


################################################################################
## Instances of family ##
################################################################################

library(taxize)

pacman::p_load(devtools, 
               tidyverse, 
               R.rsp, 
               dplyr,
               broom,
               tidyverse,
               here
)


LR.dat <- LR.dat %>%
  mutate(Study_Species = recode(Study_Species,
                                "Alces alces shirasi" = "Alces alces",
                                "Canis lupis" = "Canis lupus",
                                "Sylvilagus bachmani riparius" = "Sylvilagus bachmani",
                                "Sciurus nayaritensis chiricahuae" = "Sciurus nayaritensis",
                                "Alcelaphus buselaphus camaa" = "Alcelaphus buselaphus",
                                "Ursus arctos horribilus" = "Ursus arctos horribilis",
                                "Cyrtonix montezumae" = "Cyrtonyx montezumae"))

# Extracting Species

species_list = LR.dat$Study_Species

# Getting Class and Order

class_df <- tax_name(species_list, get = "class", db = "ncbi")
order_df <- tax_name(species_list, get = "order", db = "ncbi")

# Renaming for Join

order_df = rename(order_df, Study_Species = Species_tree)
class_df = rename(class_df, Study_Species = query)

# Making sure class_df has only one row per species, to avoid unwanted duplication

class_df_unique <- class_df %>%
  distinct(Study_Species, class)

# Join

LR.dat <- LR.dat %>%
  left_join(class_df_unique, by = "Study_Species")

# Making sure order_df has only one row per species, to avoid unwanted duplication

order_df_unique <- order_df %>%
  distinct(Study_Species, order)

# Join

LR.dat <- LR.dat %>%
  left_join(order_df_unique, by = "Study_Species")

# Naming columns to include database source

dat = rename(LR.dat, class_ncbi = class)
dat = rename(LR.dat, order_ncbi = order)

## Reducing to species level

unique(LR.dat$Study_Species)
table(LR.dat$Study_Species)

###############################################################################
## Summary stats ##

## For Each Record: Country and Year

study_metadata <- LR.dat %>%
  distinct(Firstauthor_Year, Publication_Year, Study_Country)

# Studies by year

yearly_counts <- study_metadata %>%
  count(Publication_Year, name = "n_studies") %>%
  arrange(Publication_Year)

ggplot(yearly_counts, aes(x = Publication_Year, y = n_studies)) +
  geom_col(fill = "#003366") +
  labs(
       x = "Publication Year",
       y = "Number of Studies") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),         # remove all grid lines
    axis.line = element_line(color = "black"),  # add black axis lines
    axis.ticks = element_line(color = "black"), # optional: add ticks
    axis.text = element_text(color = "black")
  )

ggsave(here("Outputs", "Figures", "Study_counts_year.pdf"))

studies_by_decade <- study_metadata %>%
  filter(!is.na(Publication_Year)) %>%
  mutate(
    Decade = floor(Publication_Year / 10) * 10,
    Decade_Label = paste0(Decade, "s")
  ) %>%
  count(Decade_Label, sort = FALSE) %>%
  mutate(
    Percent = round(n / sum(n) * 100, 1)
  ) %>%
  arrange(Decade_Label)

print(studies_by_decade)


# Cumulative 

yearly_study_counts <- study_metadata %>%
  count(Publication_Year, name = "n_studies") %>%
  arrange(Publication_Year) %>%
  mutate(cumulative_studies = cumsum(n_studies))

ggplot(yearly_study_counts, aes(x = Publication_Year, y = cumulative_studies)) +
  geom_line(color = "#2C77B2", size = 1.2) +
  geom_point(color = "#2C77B2", size = 2) +
  labs(title = "Cumulative Number of Studies Over Time",
       x = "Publication Year",
       y = "Cumulative Studies") +
  theme_minimal()

cumulative_by_country <- study_metadata %>%
  group_by(Study_Country, Publication_Year) %>%
  summarise(n_studies = n(), .groups = "drop") %>%
  arrange(Study_Country, Publication_Year) %>%
  group_by(Study_Country) %>%
  mutate(cumulative_studies = cumsum(n_studies)) %>%
  ungroup()

ggplot(cumulative_by_country, aes(x = Publication_Year, y = cumulative_studies, fill = Study_Country)) +
  geom_area(alpha = 0.8) +
  labs(title = "Cumulative Number of Studies by Country",
       x = "Publication Year",
       y = "Cumulative Studies",
       fill = "Study_Country") +
  theme_minimal()


# Studies by Country

Country_counts <- study_metadata %>%
  count(Study_Country, name = "n_studies") %>%
  arrange(Study_Country)

Country_counts_other = read.csv(here("Data/Literature Review/Country_counts_other.csv"))

# Reorder factor levels by descending count
Country_counts_other <- Country_counts_other %>%
  mutate(Study_Country = reorder(Study_Country, -n_studies))

## PLot with different colours

ggplot(Country_counts_other, aes(x = Study_Country, y = n_studies, fill = Study_Country)) +
  geom_col() +
  scale_fill_viridis_d() +
  labs(
    x = "Country",
    y = "Number of Studies"
  ) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(here("Outputs", "Figures", "Country_Count_Records.pdf"))


##plot with same colours

ggplot(Country_counts_other, aes(x = Study_Country, y = n_studies)) +
  geom_col(fill = "#003366") +
  labs(
    x = "Country",
    y = "Number of Studies") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),         # remove all grid lines
    axis.line = element_line(color = "black"),  # add black axis lines
    axis.ticks = element_line(color = "black"), # optional: add ticks
    axis.text = element_text(color = "black")
  )

ggsave(here("Outputs", "Figures", "Country_Count_Records.pdf"))

## plots by continent

continent_lookup <- tibble::tibble(
  Study_Country = c(
    "Australia", "Canada", "Indonesia", "South Africa", "South Korea",
    "Spain", "United States of America", "Kenya", "Nepal", "Paraguay",
    "Portugal", "Uganda"
  ),
  Continent = c(
    "Oceania", "North America", "Asia", "Africa", "Asia",
    "Europe", "North America", "Africa", "Asia", "South America",
    "Europe", "Africa"
  )
)

#Joining continent to country

study_metadata <- study_metadata %>%
  left_join(continent_lookup, by = "Study_Country")

# joining other study variables to create summary table for supp material

study_metadata <- study_metadata %>%
  left_join(
    LR.dat %>%
      dplyr::select(Firstauthor_Year, title, article.id, Study_Species),
    by = "Firstauthor_Year"
  ) %>%
  distinct(Firstauthor_Year, .keep_all = TRUE)



###############################################################################
###############################################################################


# Step 1: count effect sizes per Firstauthor_Year
es_counts <- es %>%
  count(Firstauthor_Year, name = "n_effect_sizes")

# Step 2: join counts into study_metadata
study_metadata <- study_metadata %>%
  left_join(es_counts, by = "Firstauthor_Year")

write.csv(study_metadata, "Outputs/supp_mat_study_metadata.csv")

# Counting studies by year and continent

study_counts <- study_metadata %>%
  group_by(Publication_Year, Continent) %>%
  summarise(n = n(), .groups = "drop")


library(viridis)

# Summarise number of studies by continent
continent_counts <- study_metadata %>%
  count(Continent) %>%
  filter(!is.na(Continent)) %>%
  mutate(Continent = fct_reorder(Continent, n, .desc = TRUE))

continent_counts <- study_metadata %>%
  count(Continent) %>%
  filter(!is.na(Continent)) %>%
  mutate(
    Percent = round(n / sum(n) * 100, 1),  # calculate percent
    Continent = fct_reorder(Continent, n, .desc = TRUE)
  )

print(continent_counts)


country_percent_within_continent <- study_metadata %>%
  filter(!is.na(Continent)) %>%
  count(Continent, Study_Country) %>%
  group_by(Continent) %>%
  mutate(
    Continent_Total = sum(n),
    Percent_of_Continent = round(n / Continent_Total * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(Continent, desc(Percent_of_Continent))

print(country_percent_within_continent)

# Plot
ggplot(continent_counts, aes(x = Continent, y = n, fill = Continent)) +
  geom_col() +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # No gap at y-axis origin
  labs(
    x = "Continent",
    y = "Number of Studies"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none"  # No legend since color = continent already shown on x
  )

ggsave(here("Outputs", "Figures", "Continent_Count_Records.pdf"))


# Multi species or impact fire in study?

multi_firetype_authors <- LR.dat %>%
  group_by(Firstauthor_Year) %>%
  summarise(n_fire_types = n_distinct(Impact_fire_type),
            fire_types = paste(sort(unique(Impact_fire_type)), collapse = ", "),
            .groups = "drop") %>%
  filter(n_fire_types > 1)

multi_firetype_authors

table(multi_firetype_authors$Firstauthor_Year)

multi_species_authors <- LR.dat %>%
  group_by(Firstauthor_Year) %>%
  summarise(n_species = n_distinct(Study_Species),
            species = paste(sort(unique(Study_Species)), collapse = ", "),
            .groups = "drop") %>%
  filter(n_species > 1)

multi_species_authors

## Summarising Study Design, Movemnet Category and Fire Type


# Study Design summary
study_design_summary <- LR.dat %>%
  count(Study_Design) %>%
  mutate(Percent = round(n / sum(n) * 100, 1))

print(study_design_summary)


LR.data_clean.mc <- LR.dat %>%
  mutate(Movement_Category_Clean = case_when(
    str_detect(Movement_category, regex("home range", ignore_case = TRUE)) &
      !str_detect(Movement_category, "&|,") ~ "Home range",
    
    str_detect(Movement_category, regex("linear distance", ignore_case = TRUE)) &
      !str_detect(Movement_category, "&|,") ~ "Linear Distance",
    
    str_detect(Movement_category, regex("tortuos", ignore_case = TRUE)) &
      !str_detect(Movement_category, "&|,") ~ "Tortuosity",
    
    str_detect(Movement_category, regex("other space", ignore_case = TRUE)) &
      !str_detect(Movement_category, "&|,") ~ "Other Space Use",
    
    str_detect(Movement_category, "&|,") ~ "Multiple",  # combos
    
    TRUE ~ Movement_category  # fallback (keep as is)
  ))

# Movement Category summary
movement_category_summary <- LR.data_clean.mc %>%
  count(Movement_Category_Clean) %>%
  mutate(Percent = round(n / sum(n) * 100, 1))

print(movement_category_summary)

# Impact Fire Type summary
impact_fire_type_summary <- LR.dat %>%
  count(Impact_fire_type) %>%
  mutate(Percent = round(n / sum(n) * 100, 1))

print(impact_fire_type_summary)


##############################################################################
## Counts by biome and fire type ######
##############################################################################

## Cleaning data, and separating lat and long into two columns

LR.dat_clean <- LR.dat %>%
  mutate(Lat_Long = str_replace_all(Lat_Long, "-\\s+", "-")) %>%
  mutate(Lat_Long = str_trim(Lat_Long)) %>%
  separate(Lat_Long, into = c("lat", "lon"), sep = ",\\s*", convert = TRUE)

# Checking for NAs

any(is.na(LR.dat_clean$lat))  
any(is.na(LR.dat_clean$lon))

## Converting to points

LR.coords_sf <- st_as_sf(LR.dat_clean, coords = c("lon", "lat"), crs = 4326)

## Pulling in Biome data

data(biomes)

## Merging biome number to dataset, using the points

LR.coords_with_biome <- st_join(LR.coords_sf, biomes, left = TRUE)

summary(LR.coords_with_biome$BIOME)

### Fixing NAs
LR.coords_with_biome %>%
  filter(is.na(BIOME))

LR.coords_with_biome %>%
  filter(is.na(BIOME)) %>%
  slice(11)

LR.coords_with_biome %>%
  filter(is.na(BIOME)) %>%
  st_coordinates()


# Plot of nas

ggplot() +
  geom_sf(data = biomes, fill = "lightgreen", color = "darkgreen") +
  geom_sf(data = LR.coords_with_biome, aes(color = is.na(BIOME))) +
  geom_text(
    data = filter(LR.coords_with_biome, is.na(BIOME)),
    aes(label = Firstauthor_Year, geometry = geometry),
    stat = "sf_coordinates",
    size = 3,
    color = "red"
  ) +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme_minimal()

LR.coords_with_biome %>%
  filter(is.na(BIOME)) %>%
  pull(Firstauthor_Year)

#############################################################################
## Creating dataframe with biome name and number

unique(LR.coords_with_biome$BIOME)

biome_no <- setdiff(1:13, 2)
biomes.names <- data.frame(
  biome_no = biome_no,
  biome_name = c(
    "Tropical and subtropical moist broadleaf forests",
    "Tropical and subtropical coniferous forests",
    "Temperate broadleaf and mixed forests",
    "Temperate coniferous forests",
    "Boreal forests/taiga",
    "Tropical and subtropical grasslands, savannas, and shrublands",
    "Temperate grasslands, savannas, and shrublands",
    "Flooded grasslands and savannas",
    "Montane grasslands and shrublands",
    "Tundra",
    "Mediterranean forests, woodlands, and scrub or sclerophyll forests",
    "Deserts and xeric shrublands"),
  biome_characteristics = c(
    "(tropical and subtropical, humid)",
    "(tropical and subtropical, semihumid)",
    "(temperate, humid)",
    "(temperate, humid to semihumid)",
    "(subarctic, humid)",
    "(tropical and subtropical, semiarid)",
    "(temperate, semiarid)",
    "(temperate to tropical, fresh or brackish water inundated)",
    "(alpine or montane climate)",
    "(Arctic)",
    "(temperate warm, semihumid to semiarid with winter rainfall)",
    "(temperate to tropical, arid)"))

biomes.names$biome_no <- as.numeric(biomes.names$biome_no)  

head(biomes.names)

LR.dat.biome <- LR.coords_with_biome %>%
  rename(biome_no = BIOME)

## Joining Biome name with number

LR.dat.biome <- LR.coords_with_biome %>%
  rename(biome_no = BIOME) %>%
  left_join(dplyr::select(biomes.names, biome_name, biome_no), by = "biome_no")


unique(LR.dat.biome$Impact_fire_type)

table(LR.dat.biome$Impact_fire_type)

count_table <- table(
  LR.dat.biome$biome_name,
  LR.dat.biome$Impact_fire_type)

count_table

# Calculate counts of fire types per biome
plot_data <- LR.dat.biome %>%
  count(biome_name, Impact_fire_type, name = "Count")


biome_totals <- plot_data %>%
  group_by(biome_name) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Impact_fire_type = "Total")

# Combine totals with original data
plot_data_total <- bind_rows(plot_data, biome_totals)

print(plot_data_total)

# Convert to factors, include "Total" as last level for Impact_fire_type
plot_data_total <- plot_data_total %>%
  mutate(
    biome_name = factor(biome_name),
    Impact_fire_type = factor(Impact_fire_type, levels = c(levels(factor(plot_data$Impact_fire_type)), "Total"))
  )

ggplot(plot_data_total, 
       aes(x = as.numeric(Impact_fire_type), 
           y = as.numeric(biome_name), 
           fill = Count)) +
  
  geom_tile(width = 1, height = 1, color = NA) +
  
  geom_text(aes(label = Count,
                color = ifelse(Count > max(Count)/2, "white", "black")),
            size = 3.5) +
  scale_color_identity() +
  
  scale_fill_gradient(low = "#E6F2FF", high = "#003366", name = "Instances") +
  
  scale_x_continuous(
    breaks = seq_along(levels(plot_data_total$Impact_fire_type)),
    labels = levels(plot_data_total$Impact_fire_type),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(plot_data_total$biome_name)),
    labels = levels(plot_data_total$biome_name),
    expand = c(0, 0)
  ) +
  
  labs(x = "Fire Type", y = "Biome") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks.length = unit(5, "pt"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.line = element_line(color = "black", size = 0.7),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# Saving Plot

ggsave(here("Outputs", "Figures", "LR_CountbyBiome.pdf"))


biome_summary <- LR.dat.biome %>%
  count(biome_name, name = "Count") %>%
  mutate(
    Percent = round(Count / sum(Count) * 100, 1)
  ) %>%
  arrange(desc(Count))

print(biome_summary)

################################################################################
## Instances of family ##
################################################################################

library(taxize)

pacman::p_load(devtools, 
               tidyverse, 
               R.rsp, 
               dplyr,
               broom,
               tidyverse,
               here
)


LR.dat <- LR.dat %>%
  mutate(Study_Species = recode(Study_Species,
                                "Alces alces shirasi" = "Alces alces",
                                "Canis lupis" = "Canis lupus",
                                "Sylvilagus bachmani riparius" = "Sylvilagus bachmani",
                                "Sciurus nayaritensis chiricahuae" = "Sciurus nayaritensis",
                                "Alcelaphus buselaphus camaa" = "Alcelaphus buselaphus",
                                "Ursus arctos horribilus" = "Ursus arctos horribilis",
                                "Cyrtonix montezumae" = "Cyrtonyx montezumae"))

# Extracting Species

species_list = LR.dat$Study_Species

# Getting Class and Order

class_df <- tax_name(species_list, get = "class", db = "ncbi")
order_df <- tax_name(species_list, get = "order", db = "ncbi")
family_df = tax_name(species_list, get = "family", db = "ncbi")

# Renaming for Join

order_df = rename(order_df, Study_Species = query)
class_df = rename(class_df, Study_Species = query)
family_df = rename(family_df, Study_Species = query)

# Making sure class_df has only one row per species, to avoid unwanted duplication

class_df_unique <- class_df %>%
  distinct(Study_Species, class)

# Join

LR.dat <- LR.dat %>%
  left_join(class_df_unique, by = "Study_Species")

# Making sure order_df has only one row per species, to avoid unwanted duplication

order_df_unique <- order_df %>%
  distinct(Study_Species, order)

# Join

LR.dat <- LR.dat %>%
  left_join(order_df_unique, by = "Study_Species")

# Making sure familyr_df has only one row per species, to avoid unwanted duplication

family_df_unique <- family_df %>%
  distinct(Study_Species, family)

# Join

LR.dat <- LR.dat %>%
  left_join(family_df_unique, by = "Study_Species")

table(LR.dat$family)

table(LR.dat$order)

missing_class <- LR.dat %>%
  filter(is.na(class)) %>%
  distinct(order)

print(missing_class)
LR.dat$class[LR.dat$order == "Testudines" & is.na(LR.dat$class)] <- "Reptilia"
LR.dat$class[LR.dat$order == "Squamata"] <- "Reptilia"

LR.dat$class <- as.factor(LR.dat$class)
LR.dat$order <- as.factor(LR.dat$order)

LR.dat <- LR.dat %>%
  mutate(order_grouped = fct_reorder(order, as.numeric(factor(class))))

ggplot(LR.dat, aes(x = order_grouped, fill = class)) +
  geom_bar() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(option = "D") +
  labs(
    x = "Order",
    y = "Count",
    fill = "Class"
  )

ggsave(here("Outputs", "Figures", "LR_CountbyOrder.pdf"))

order_class_summary <- LR.dat %>%
  count(class, order, name = "Count") %>%
  group_by(class) %>%
  mutate(Percent_within_class = round(100 * Count / sum(Count), 1)) %>%
  ungroup() %>%
  arrange(class, desc(Percent_within_class))

print(order_class_summary)

order_summary <- LR.dat %>%
  count(order, name = "Count") %>%
  mutate(Percent = round(100 * Count / sum(Count), 1)) %>%
  arrange(desc(Percent))

print(order_summary)

class_summary <- LR.dat %>%
  count(class, name = "Count") %>%
  mutate(Percent = round(100 * Count / sum(Count), 1)) %>%
  arrange(desc(Percent))

print(class_summary)

species_summary = LR.dat %>%
  count(Study_Species, name = "Count") %>%
  mutate(Percent = round(100 * Count / sum(Count), 1 )) %>%
  arrange(desc(Percent))

print(species_summary)      
