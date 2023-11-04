
# Libraries

library(tidyverse)
library(sf)


# Data

grid <- st_read("data/spatial/mod/boundaries/spatial_block_grid_w_attr.gpkg") %>%
  rename(NoSnowDOY = calDoy_median,
         Greenup_Date = Greenup_1,
         MidGreenup_Date = MidGreenup_1,
         Maturity_Date = Maturity_1,
         Peak_Date = Peak_1,
         MidGreendown_Date = MidGreendown_1,
         Senescence_Date = Senescence_1,
         Dormancy_Date = Dormancy_1) %>%
  select(-c(relDoy_median,count,label,id)) %>%
  mutate_at(vars(contains('_Date')), 
            function(x) (as.Date(as.POSIXct(x*24*60*60, origin = "1970-01-01", tz="UTC")))) %>%
  mutate(group = "1")

glimpse(grid)

# Grab some statistics across all grids

median(grid$NoSnowDOY)
median(as.numeric(strftime(grid$Greenup_Date, format = "%j")))
median(as.numeric(strftime(grid$MidGreenup_Date, format = "%j")))
median(as.numeric(strftime(grid$Maturity_Date, format = "%j")))
median(as.numeric(strftime(grid$MidGreendown_Date, format = "%j")))
median(as.numeric(strftime(grid$Senescence_Date, format = "%j")))
median(as.numeric(strftime(grid$Dormancy_Date, format = "%j")))


# Boxplot

# Reshape
grid.m <- reshape2::melt(
  grid,id.vars='group', 
  measure.vars=c(
    'Greenup_Date','MidGreenup_Date','Maturity_Date',
    'Peak_Date','MidGreendown_Date','Senescence_Date','Dormancy_Date'))

p <- ggplot(grid.m) +
  geom_boxplot(aes(x=group, y=value, color=variable)) +
  scale_y_date(date_labels="%b",date_breaks  ="1 month") +
  labs(x="",y="Date") +
  coord_flip() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
p

# Spatial maps
ggplot() +
  geom_sf(data=grid, aes(fill=Greenup_Date)) +
  theme_void()

# Export the shapefile
st_write(grid,"data/spatial/mod/boundaries/spatial_cv_grid_1km_w_attr.shp",delete_dsn=TRUE)

