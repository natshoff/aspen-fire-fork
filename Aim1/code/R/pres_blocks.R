
library(tidyverse)
library(sf)

# Load the grids, keep the elevation attribute and TreeMap summary
grid <- st_read('data/spatial/mod/boundaries/spatial_block_grid_50km2.gpkg') %>%
  select(grid_id,elevation_mn,treemap_sum) %>%
  mutate(treemap_m2 = treemap_sum*900,
         grid_area = as.numeric(st_area(geom)),
         treemap_pct = (treemap_m2 / grid_area)*100)
glimpse(grid)

# Spatial map of aspen %
ggplot() +
  geom_sf(data=grid, aes(fill=treemap_pct)) +
  scale_fill_viridis_c(option="magma") +
  theme_light()

# Get a summary of aspen area percent across grids
summary(grid$treemap_pct)

# How many grids have at least 5% aspen area
dim(grid %>% filter(treemap_pct >= 5))


############################################

# Load the photo interpretation points
Pres <- st_read('data/spatial/mod/training/points/gee/pi_points_srme_m500_.shp') %>%
  st_transform(st_crs(grid)) %>%
  # Join to the spatial block grid
  st_join(grid)
glimpse(Pres)

# Get a count of presence data by grid id
count <- Pres %>%
  filter(treemap_pct >= 5) %>%
  group_by(grid_id) %>%
  summarize(pres_count = n())

# Plot the histogram
ggplot(data=count, aes(x=pres_count)) +
  geom_histogram()

summary(count$pres_count)

# How many have less than 100 points?
dim(count %>% filter(pres_count >= 100))
dim(count %>% filter(pres_count >= 50))

# count_l100 <- count %>% filter(pres_count < 100)
# # Identify blocks which have less than 100 and export as a feature
# grid_l100 <- grid %>% filter(grid_id %in% count_l100$grid_id) %>%
#   # Join the presence count
#   left_join(count_l100%>%select(grid_id,pres_count)%>%as_tibble(), by="grid_id")
# st_write(grid_l100,'data/spatial/mod/training/spatial_block_grid_50km_less100addin.gpkg')
