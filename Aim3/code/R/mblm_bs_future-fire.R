# Calculate Thiel-Sen estimators for temporal trends in future fire
# By US Ecoregion Level IV (from "Fires of Unusual Size")

# libraries
library(tidyverse)
library(sf)
library(mblm)
library(viridis)
library(patchwork)
library(ggspatial)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load the data (future fire)
ff <- read.csv(paste0(projdir,'data/tabular/us_l4eco_future-fire.csv'))
glimpse(ff)

# make sure the data frame is distinct
ff <- ff %>% distinct(US_L4CODE, Year, .keep_all = TRUE)

# Function to apply MBLM (Thiel-Sen estimator) per ecoregion
ts_trends <- function(df) {
 # Fit Theil-Sen models
 area_trend <- summary(mblm(En_Area ~ Year, data = df, repeated = F))
 count_trend <- summary(mblm(En_NFire ~ Year, data = df, repeated = F))
 # Extract slopes (trend per year)
 data.frame(
  US_L4CODE = df$US_L4CODE[1],
  trend_area = area_trend$coefficients[2, 1],  # Theil-Sen slope for fire size
  p_area = area_trend$coefficients[2, 4],      # P-value for fire size trend
  trend_count = count_trend$coefficients[2, 1], # Theil-Sen slope for fire count
  p_count = count_trend$coefficients[2, 4]      # P-value for fire count trend
 )
}

# run the mblm for ecoregions
eco.trends <- ff %>%
 group_by(US_L4CODE) %>%
 group_split() %>%
 map_dfr(ts_trends)

# Display results
print(eco.trends)


##################################################
# load the spatial data and merge the results back
ecol4 <- st_read(paste0(projdir,'data/spatial/raw/boundaries/srm_eco_l4.gpkg'))
glimpse(ecol4)

# join the results
ecol4 <- ecol4 %>%
 left_join(eco.trends, by='US_L4CODE')
glimpse(ecol4)


#################
# fire area trend
p1 <- ggplot(ecol4) +
 geom_sf(aes(fill = trend_area), color = NA) +  # No border for smooth visualization
 scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Burned area trend") +
 # scale_fill_viridis(option = "inferno", direction = -1, name = "Fire Size Trend") +
 guides(fill = guide_colourbar(direction = "vertical", 
                               barwidth = 0.8, 
                               barheight = 10, 
                               ticks.colour = NA, 
                               title.position = "left")) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.25, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.15,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.25,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p1

# fire counts trend
p2 <- ggplot(ecol4) +
 geom_sf(aes(fill = trend_count), color = NA) +  # No border for smooth visualization
 scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Fire occurrence trend") +
 # scale_fill_gradient(low = "#ffedea", high = "#b30000") + 
 # scale_fill_viridis(option = "YlOrRd", direction = 1, name = "Fire Count Trend") +
 guides(fill = guide_colourbar(direction = "vertical", 
                               barwidth = 0.8, 
                               barheight = 10, 
                               ticks.colour = NA, 
                               title.position = "left")) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.25, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.15,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.25,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)

# combine the two plots
p3 <- p1 + p2
p3

# save the plot.
out_png <- paste0(projdir,'figures/SouthernRockies_FutureFire_Trends.png')
ggsave(out_png, plot = p3, dpi = 500, height=6.5, width=7, bg = 'white')

# save the spatial data
out_fp <- paste0(projdir,'data/spatial/mod/future_fire_ecol4_trend.gpkg')
st_write(ecol4, out_fp, append = F)



