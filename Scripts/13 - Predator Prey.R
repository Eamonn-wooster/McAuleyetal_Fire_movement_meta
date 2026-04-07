

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
               sf,
               geodat,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp #,#added
               # miWQS #added
) 



# 1. Load the datasets from Hao et al.
metaweb <- readRDS(here("Data/Traits/predicted.links.matrix.RDS"))              # global predator–prey pairs
grid_metrics <- readRDS(here("Data/Traits/Tetrapods_Network_Metrics_with_coordinate.RDS")) # grid IDs + coords

dat = read.csv(here("Outputs/old/Dat_for_analysis-ss-aves-rep-mam-mass-Rev3.csv"))

dat_clean <- dat %>%
  mutate(Lat_Long = str_replace_all(Lat_Long, "-\\s+", "-")) %>%
  mutate(Lat_Long = str_trim(Lat_Long)) %>%
  separate(Lat_Long, into = c("lat", "lon"), sep = ",\\s*", convert = TRUE)

# Checking for NAs

any(is.na(dat_clean$lat))  
any(is.na(dat_clean$lon))

## Converting to points

coords_sf <- st_as_sf(dat_clean, coords = c("lon", "lat"), crs = 4326)


# 3. Match study location to nearest Hao et al. grid cell
grid_sf <- st_as_sf(grid_metrics, coords = c("x", "y"), crs = 4326)
study_sf <- st_as_sf(coords_sf, coords = c("lon", "lat"), crs = 4326)
study_sf <- st_join(coords_sf, grid_sf["grid_id"], join = st_nearest_feature)

any(is.na(study_sf$grid_id))  

# 4. Subset the metaweb to species present in that grid cell
# (if you have a list of species for that study, replace below with your actual list)
##local_species <- c("Macropus rufus", "Panthera leo", "Canis lupus", "Rattus rattus")

local_species <- data.frame(study_sf$Species)

local_links <- metaweb %>%
  filter(predator %in% local_species & prey %in% local_species)

# 5. Build a local network and compute prey-averaged trophic levels (PATL)
g <- graph_from_data_frame(local_links, directed = TRUE)

# Function for recursive trophic-level computation
calc_PATL <- function(g) {
  basal <- V(g)[degree(g, mode = "out") == 0]
  trophic <- rep(NA, vcount(g))
  names(trophic) <- V(g)$name
  trophic[basal$name] <- 1
  repeat {
    uncalculated <- is.na(trophic)
    if (!any(uncalculated)) break
    for (v in names(trophic[uncalculated])) {
      prey <- neighbors(g, v, mode = "out")$name
      if (length(prey) > 0 && all(!is.na(trophic[prey]))) {
        trophic[v] <- 1 + mean(trophic[prey])
      }
    }
  }
  return(trophic)
}

PATL <- calc_PATL(g)
patl_df <- data.frame(species = names(PATL), PATL = PATL)

# 6. Save or merge with your meta-analysis dataset
write.csv(patl_df, "local_PATL_example.csv", row.names = FALSE)

# You can now merge this PATL dataframe with your movement dataset
# and include PATL (trophic position) as a moderator in your models.

