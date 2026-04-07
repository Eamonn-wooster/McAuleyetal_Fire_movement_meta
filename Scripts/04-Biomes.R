## Using geodat package (https://github.com/RS-eco/geodat/blob/main/R/geodat-package.R) ##

#' @docType data
#' @name biomes
#' @title Biomes of the world
#' @description Shapefile of the world's biomes extracted from Olson's ecoregions of the worlds.
#' @usage data(biomes)
#' @details This shapefile contains information about the world's biomes and 
#' was derived from the terrestrial ecoregions of the world'a shapefile.
#' 

# Packages #
remotes::install_github("RS-eco/geodat", force = TRUE)
library(geodat)
library(sf)     
library(dplyr) 
library(tidyr)    
library(stringr)  
library(terra) 
library(here)
library(readr)

## Importing extraction data

dat = read.csv(here("Outputs/Dat_for_analysis-ss-mass-Rev4.csv"))

## Cleaning data, and separating lat and long into two columns

dat_clean <- dat %>%
  mutate(Lat_Long = str_replace_all(Lat_Long, "-\\s+", "-")) %>%
  mutate(Lat_Long = str_trim(Lat_Long)) %>%
  separate(Lat_Long, into = c("lat", "lon"), sep = ",\\s*", convert = TRUE)

# Checking for NAs

any(is.na(dat_clean$lat))  
any(is.na(dat_clean$lon))

## Converting to points

coords_sf <- st_as_sf(dat_clean, coords = c("lon", "lat"), crs = 4326)

## Pulling in Biome data

data(biomes)

## Merging biome number to dataset, using the points

coords_with_biome <- st_join(coords_sf, biomes, left = TRUE)

summary(coords_with_biome$BIOME)

### Fixing NAs

# 1 Identify points with NA biome
problem_points <- coords_with_biome %>%
  filter(is.na(BIOME))

# 2 Find nearest biome for each NA point
nearest_index <- st_nearest_feature(problem_points, biomes)

# 3 Extract the nearest biome ID
nearest_biomes <- biomes[nearest_index, ]$BIOME

# 4 Assign to the NA points in the main dataset
coords_with_biome$BIOME[is.na(coords_with_biome$BIOME)] <- nearest_biomes

summary(coords_with_biome$BIOME)

## Creating dataframe with biome name and number

biomes.names <- data.frame(
  biome_no = 1:14,
  biome_name = c(
    "Tropical and subtropical moist broadleaf forests",
    "Tropical and subtropical dry broadleaf forests",
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
    "Deserts and xeric shrublands",
    "Mangrove"),
  biome_characteristics = c(
  "(tropical and subtropical, humid)",
  "(tropical and subtropical, semihumid)",
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
  "(temperate to tropical, arid)",
  "(subtropical and tropical, salt water inundated)"))

biomes.names$biome_no <- as.numeric(biomes.names$biome_no)  

head(biomes.names)

dat.biome = coords_with_biome %>%
  dplyr::rename(biome_no = BIOME)

## Joining Biome name with number

dat.biome <- dat.biome %>%
  left_join(biomes.names %>% dplyr::select(biome_name, biome_no),
            by = "biome_no")

## Writing CSV

write_csv(dat.biome, "Outputs/Dat_for_analysis-ss-mass-biome-Rev4.csv")


################################## Biome to veg_climate #######################

dat.biome <- read.csv(here::here("Outputs/Dat_for_analysis-ss-mass-biome-Rev4.csv")) 

## Vegetation structure

table(dat.biome$biome_name)

dat.biome$veg_str <- dplyr::case_when(
  dat.biome$biome_name %in% c(
    "Boreal forests/taiga",
    "Mediterranean forests, woodlands, and scrub or sclerophyll forests",
    "Temperate broadleaf and mixed forests",
    "Temperate coniferous forests",
    "Tropical and subtropical coniferous forests",
    "Tropical and subtropical moist broadleaf forests"
  ) ~ "closed",
  
  dat.biome$biome_name %in% c(
    "Deserts and xeric shrublands",
    "Flooded grasslands and savannas",
    "Montane grasslands and shrublands",
    "Temperate grasslands, savannas, and shrublands",
    "Tropical and subtropical grasslands, savannas, and shrublands"
  ) ~ "open",
  
  TRUE ~ NA_character_
)

dat.biome$veg_str <- factor(dat.biome$veg_str, levels = c("closed", "open"))

table(dat.biome$veg_str)

any(is.na(dat.biome$veg_str)) 

dat.biome %>%
  filter(is.na(veg_str))

## Climate

## Cleaning data, and separating lat and long into two columns

dat.biome$geometry <- gsub("c\\(|\\)", "", dat.biome$geometry)

dat_clean.clim <- dat.biome %>%
  mutate(geometry = str_replace_all(geometry, "-\\s+", "-")) %>%
  mutate(geometry = str_trim(geometry)) %>%
  separate(geometry, into = c("lon", "lat"), sep = ",\\s*", convert = TRUE)

# Checking for NAs

any(is.na(dat_clean.clim$lat))  
any(is.na(dat_clean.clim$lon))

# to vector 

pts_vect <- vect(dat_clean.clim)

# Raster from Beck 2023

kg_raster <- rast(here::here("Data/Climate/1991_2020/1991_2020/koppen_geiger_0p1.tif"))

# Check coordinate reference system (CRS)
crs(kg_raster)

dat_clean.clim$climate_zone <- terra::extract(kg_raster, pts_vect)[, 2]

any(is.na(dat_clean.clim$climate_zone)) 

# Adding zone names

koppen_lookup <- data.frame(
  climate_zone = 1:30,
  code = c(
    "Af","Am","Aw","BWh","BWk","BSh","BSk",
    "Csa","Csb","Csc","Cwa","Cwb","Cwc",
    "Cfa","Cfb","Cfc","Dsa","Dsb","Dsc","Dsd",
    "Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd",
    "ET","EF"
  ),
  name = c(
    "Tropical, rainforest",
    "Tropical, monsoon",
    "Tropical, savannah",
    "Arid, desert, hot",
    "Arid, desert, cold",
    "Arid, steppe, hot",
    "Arid, steppe, cold",
    "Temperate, dry summer, hot summer",
    "Temperate, dry summer, warm summer",
    "Temperate, dry summer, cold summer",
    "Temperate, dry winter, hot summer",
    "Temperate, dry winter, warm summer",
    "Temperate, dry winter, cold summer",
    "Temperate, no dry season, hot summer",
    "Temperate, no dry season, warm summer",
    "Temperate, no dry season, cold summer",
    "Cold, dry summer, hot summer",
    "Cold, dry summer, warm summer",
    "Cold, dry summer, cold summer",
    "Cold, dry summer, very cold winter",
    "Cold, dry winter, hot summer",
    "Cold, dry winter, warm summer",
    "Cold, dry winter, cold summer",
    "Cold, dry winter, very cold winter",
    "Cold, no dry season, hot summer",
    "Cold, no dry season, warm summer",
    "Cold, no dry season, cold summer",
    "Cold, no dry season, very cold winter",
    "Polar, tundra",
    "Polar, frost"
  ),
  stringsAsFactors = FALSE
)

# Joining 

dat_clean.clim <- dat_clean.clim %>%
  left_join(koppen_lookup, by = "climate_zone")

# Creating Broad Categories

dat_clean.clim <- dat_clean.clim %>%
  mutate(Climate = str_replace(word(name, 1), ",", ""))

# Creating Veg_climate category

dat_clean.clim <- dat_clean.clim %>%
  mutate(veg_climate = paste(veg_str, Climate, sep = "_"))

dat_clean.clim$veg_climate = as.factor(dat_clean.clim$veg_climate)

table(dat_clean.clim$veg_climate)

## Looking into closed_Arid studies

dat_clean.clim %>%
  filter(veg_climate == "closed_Arid") #both from high elevation in Arizona, seems legit
#although Gwinn2016 should likely be Temperate coniferous forests, not Tropical and subtropical coniferous forests.

## Looking into open_Temperate studies

dat_clean.clim %>%
  filter(veg_climate == "open_Temperate")
# 4 cases from 1 Tropical and subtropical grasslands, savannas, and shrublands...

## Writing CSV

write_csv(dat_clean.clim, "Outputs/Dat_for_analysis-ss-mass-biome-clim-Rev4.csv")

## Adding Fire Activity

# Importing fire activity shapefile

fireact <- st_read(here("Data/PausasRibeiro2017GEBmaps/Pausas-Ribeiro-2017-GEB-fire-div-maps.shp"))

coords_fireact <- st_join(coords_sf, fireact)

any(is.na(coords_fireact$fireactivi)) 

summary(coords_fireact$fireactivi)

plot(st_geometry(fireact), col = 'lightgrey', border = NA)
plot(st_geometry(coords_sf), col = 'black', add = TRUE)
plot(st_geometry(coords_fireact %>% filter(is.na(fireactivi))), col = 'red', add = TRUE)

coords_fireact <- st_join(
  st_buffer(coords_sf, dist = 150),   # buffer to capture NAs
  fireact
)

any(is.na(coords_fireact$fireactivi)) 

summary(coords_fireact$fireactivi)

## Joining to main data

dat.all <- dat_clean.clim %>%
  left_join(
    coords_fireact %>%
      dplyr::select(Study_ID, ecoreg, diversity, fireactivi),
    by = "Study_ID"
  )

write_csv(dat.all, "Outputs/Dat_for_analysis-ss-mass-biome-clim-fireact-Rev4.csv")
