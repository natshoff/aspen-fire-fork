
# Load the libraries
library(tidyverse)
library(sf)
library(biscale)
library(cowplot)
library(viridis)
library(ggspatial)
library(gridExtra)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/"

# Load the data (firesheds)
firesheds <- st_read(paste0(projdir, 'Aim3/data/spatial/mod/srm_firesheds_model_data.gpkg'))
glimpse(firesheds)


######################################################
#============= aspen suitability map ================#

p1 <- ggplot(firesheds) +
 geom_sf(aes(fill = historic), color = NA) +  # No border for smooth visualization
 # scale_fill_viridis(option = "viridis", na.value = "gray80", 
 #                    name = "Historic Aspen Suitability") +
 scale_fill_distiller(palette = "Greens", direction = 1, 
                      name = "Historic aspen suitability") +
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
# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_FutureAspen_Historic.png')
ggsave(out_png, plot = p1, dpi = 500, bg = 'white')


################################
# change in future suitability #
p2 <- ggplot(firesheds) +
 geom_sf(aes(fill = delta585), color = NA) +  # No border for smooth visualization
 # scale_fill_viridis(option = "viridis", na.value = "gray80", 
 #                    name = "Change in Aspen Suitability") +
 scale_fill_distiller(palette = "RdBu", direction = 1, 
                      name = "Aspen Suitability Change") +
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
# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_FutureAspen_SSP585.png')
ggsave(out_png, plot = p2, dpi = 500, bg = 'white')

######################
# merge the two maps #
p_merged <- grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ncol = 2)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_FutureAspen_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 8, height = 6, bg = "white")


################################################
#============= future fire map ================#

p1 <- ggplot(firesheds) +
 geom_sf(aes(fill = trend_count), color = NA) +  # No border for smooth visualization
 # scale_fill_viridis(option = "viridis", na.value = "gray80", 
 #                    name = "Historic Aspen Suitability") +
 scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                      name = "Fire occurrence trend") +
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

################################
# change in future suitability #
p2 <- ggplot(firesheds) +
 geom_sf(aes(fill = trend_area), color = NA) +  # No border for smooth visualization
 # scale_fill_viridis(option = "viridis", na.value = "gray80", 
 #                    name = "Change in Aspen Suitability") +
 scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                      name = "Burned area trend") +
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

######################
# merge the two maps #
p_merged <- grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ncol = 2)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_FutureFire_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 8, height = 6, bg = "white")


######################################################
#============= built environment maps================#

####################
# WUI designations #

# calculate the sum of WUI classes
firesheds <- firesheds %>%
 mutate(wui_sum = wui1 + wui2 + wui3 + wui4)

p1 <- ggplot(firesheds) +
 geom_sf(data = firesheds %>% filter(wui_sum == 0), fill = "gray90", color = NA) + 
 geom_sf(data = firesheds %>% filter(wui_sum > 0), aes(fill = wui_sum), color = NA) +  
 scale_fill_viridis(option = "rocket", name = "Total WUI %", direction=-1) +
 # scale_fill_distiller(palette = "Greens", direction = 1, 
 #                      name = "Total WUI %", na.value = "gray80") +
 guides(fill = guide_colourbar(
  direction = "vertical", 
  barwidth = 0.4, 
  barheight = 6, 
  ticks.colour = NA, 
  title.position = "left")
 ) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p1
# # save the plot.
# out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_WUI_total.png')
# ggsave(out_png, plot = p2, dpi = 500, bg = 'white')


#################
# COMBUST (sum) #

p2 <- ggplot(firesheds) +
 geom_sf(data = firesheds %>% filter(combust_sum == 0), fill = "gray90", color = NA) + 
 geom_sf(data = firesheds %>% filter(combust_sum > 0), aes(fill = combust_sum), color = NA) + 
 scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                    name = "COMBUST (sum)", trans="sqrt",
                    labels = scales::label_number(scale = 1e-4, suffix = "k")) +
 guides(fill = guide_colourbar(
  direction = "vertical", 
  barwidth = 0.4, 
  barheight = 6, 
  ticks.colour = NA, 
  title.position = "left")
 ) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p2


###################
# Distance to WUI #

p3 <- ggplot(firesheds) +
 geom_sf(data = firesheds, aes(fill = wui_dist_mean), color = NA) + 
 scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                    name = "Distance to WUI") +
 # scale_fill_distiller(palette = "Greens", direction = 1, 
 #                      name = "Historic aspen suitability") +
 guides(fill = guide_colourbar(direction = "vertical", 
                               barwidth = 0.4, 
                               barheight = 6, 
                               ticks.colour = NA, 
                               title.position = "left")) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p3


######################
# merge the two maps #
p_merged <- grid.arrange(
 ggplotGrob(p1), 
 ggplotGrob(p2), 
 ggplotGrob(p3), 
 ncol = 3
)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_BuiltEnv_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 9, height = 6, bg = "white")


################################################
#============= vegetation maps ================#

# Reshape data to long format
firesheds_l <- firesheds %>%
 pivot_longer(cols = c(Lodgepole, Aspen, Ponderosa, Piñon_juniper, 
                       Sagebrush, Spruce_fir, Douglas_fir, White_fir, Gambel_oak),
              names_to = "forest", values_to = "pct_cover")

# create a facet map
p1 <- ggplot(firesheds_l) +
 geom_sf(aes(fill = pct_cover), color = NA) +
 scale_fill_distiller(palette = "Greens", direction = 1,
                      name = "Percent cover", na.value = "gray80") +
 facet_wrap(~forest, nrow = 2) +  # Creates small multiples by forest type
 theme_void() +
 theme(legend.position = c(0.95, 0.26),
       legend.title = element_text(angle = 90, size = 10),
       legend.key.width = unit(4, "cm"),
       legend.key.height = unit(0.4, "cm"),
       strip.text = element_text(size = 9, face = "italic",
                                 margin = margin(b = 2)),
       panel.spacing = unit(0.5, "lines")) +
 guides(fill = guide_colorbar(
  ticks.colour = NA, 
  title.position = "left",
  title.hjust = 0.5,  
  barwidth = unit(0.65, "cm"), 
  barheight = unit(6.5, "cm")  
 ))
p1

# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_EVTSAF_Panel.png")
ggsave(out_png, plot = p1, dpi = 500, width = 9, height = 6, bg = "white")


