# spatial grid
fp <- paste0(maindir,"data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_aspenfires_pixar_gridstats.gpkg")
grid <- st_read(fp)
glimpse(grid)

# Load the fire perimeter dataset
fp <- paste0(maindir,"data/spatial/mod/NIFC/nifc-ics_2018_to_2023-aspen.gpkg")
fires <- st_read(fp) %>%
 rename(Fire_ID = NIFC_ID) %>%
 filter(Fire_ID %in% grid$Fire_ID) %>%
 select(Fire_ID, geom)
glimpse(fires)

# Distance to fire center for spatial structure

fire_centers <- st_centroid(fires)

st_crs(grid) == st_crs(fire_centers)

grid <- grid %>%
 left_join(fire_centers %>% st_drop_geometry(), by = "Fire_ID") %>%
 mutate(center_geom = st_geometry(fire_centers)[match(Fire_ID, fire_centers$Fire_ID)])

# Compute distances to the fire center
grid$fire_dist <- st_distance(st_geometry(grid), st_sfc(grid$center_geom, crs = st_crs(grid)))
grid$fire_dist <- as.numeric(grid$fire_dist)  # Convert to numeric
colnames(grid)

rm(grid)
gc()