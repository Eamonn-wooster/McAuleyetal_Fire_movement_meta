library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggnewscale)

# Get world basemap
world <- ne_countries(scale = "small", returnclass = "sf")

# get data

data <- read.csv(here("Outputs/Dat_for_analysis-LR-biome.csv"))

### Extract study level metadata from case data and tidy Impact Fire Type 

# getting study level info on impact fire types for studies that looked at multiple impact fire types

data <- data %>% 
  group_by(Firstauthor_Year) %>%
  mutate(
    Impact_Fire_Types = paste(unique(Impact_fire_type), collapse = ", ")
  ) %>%
  ungroup()

study_metadata <- data%>%
  distinct(Firstauthor_Year, Publication_Year, Study_Country, Impact_Fire_Types)

study_metadata <- study_metadata %>%
  left_join(
    data %>%
      dplyr::select(Firstauthor_Year, title, article.id, Study_Species, geometry),
    by = "Firstauthor_Year"
  ) %>%
  distinct(Firstauthor_Year, .keep_all = TRUE)

study_metadata$Impact_Fire_Types

study_metadata <- study_metadata %>%
  mutate(Impact_Fire_Types = recode(Impact_Fire_Types,
                                    "Both, Planned, Wildfire" = "Both"))

# Extract long/lat from the geometry column
metadata_clean <- study_metadata %>%
  mutate(
    geometry = gsub("c\\(|\\)", "", geometry),   # remove "c(" and ")"
    lon = as.numeric(sub(",.*", "", geometry)),  # part before comma
    lat = as.numeric(sub(".*,", "", geometry))   # part after comma
  )

# Convert to sf object
metadata_sf <- st_as_sf(metadata_clean, coords = c("lon", "lat"), crs = 4326)

#Getting counts of countries

Country_counts <- study_metadata %>%
  count(Study_Country, name = "n_studies") %>%
  arrange(Study_Country)

Country_counts <- Country_counts %>%
  rename(name = Study_Country)


#joining to world base map to plot
world_counts <- left_join(world, Country_counts, by = "name") %>%
  mutate(n_studies = ifelse(is.na(n_studies), 0, n_studies))

# Converting to numeric and factor

world_counts$n_studies <- as.numeric(world_counts$n_studies)
world_counts$n_studies_fill <- ifelse(world_counts$n_studies == 0, NA, world_counts$n_studies)

metadata_sf$Impact_Fire_Types <- as.factor(metadata_sf$Impact_Fire_Types)

summary(metadata_sf$Impact_Fire_Types)

## Plot!

#mercator projection

ggplot() +
  # Countries fill with transparency
  geom_sf(data = world_counts, aes(fill = n_studies_fill),
          color = NA, alpha = 0.7) +
  # Country outlines
  geom_sf(data = world_counts, fill = NA,
          color = "black", size = 0.05) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "log",
    name = "Studies per Country",
    na.value = alpha("lightgrey", 0.2), 
    breaks = scales::log_breaks(n = 5)
  ) +
  
  # Start new fill scale for points
  new_scale_fill() +
  
  # Points layer
  geom_sf(data = metadata_sf, aes(fill = Impact_Fire_Types),
          shape = 21, color = "black", size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Both" = "#ff4d6d",
      "Unclear" = "green",
      "Planned" = "#fde725",
      "Wildfire" = "#c77cff",
      name = "Impact fire type")) +
  
  # Gray theme with legend on the right
  theme_gray() +
  theme(
    legend.position = "right",
    legend.key = element_blank() 
  )

# Robinson projection

robinson_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Transform both world and points
world_robinson <- st_transform(world_counts, crs = robinson_crs)
metadata_sf_robinson <- st_transform(metadata_sf, crs = robinson_crs)

summary_map = ggplot() +
  # Countries fill with transparency
  geom_sf(data = world_robinson, aes(fill = n_studies_fill),
          color = NA, alpha = 0.7) +
  # Country outlines
  geom_sf(data = world_robinson, fill = NA,
          color = "black", size = 0.05) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "log",
    name = "Studies per Country",
    na.value = alpha("lightgrey", 0.2), 
    breaks = scales::log_breaks(n = 5)
  ) +
  
  # Start new fill scale for points
  new_scale_fill() +
  
  # Points layer
  geom_sf(data = metadata_sf_robinson, aes(fill = Impact_Fire_Types),
          shape = 21, color = "black", size = 2, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "Both" = "#ff4d6d",
      "Unclear" = "green",
      "Planned" = "#fde725",
      "Wildfire" = "#c77cff"), 
      name = "Impact Fire Type") +
  
  # Gray theme with legend on the right
  theme_void() +
  theme(
    legend.position = "right",
    legend.key = element_blank() 
  )

summary_map ## some countries are filled with no corresponding points, need to sort out

ggsave(here("Outputs/Figures", "summary map.pdf"), height = 4, width = 8)



