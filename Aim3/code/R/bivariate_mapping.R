
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



############################################
# Bivariate map: Future Fire -> Future Aspen

# scale the variables
firesheds.biv <- firesheds %>%
 mutate(across(c(trend_count, delta585), ~scale(.)))

# create the bivariate classes
firesheds.biv <- bi_class(
 firesheds, x = trend_count, y = delta585, 
 style = "quantile", dim = 3)

# create the bivariate chloropleth map
p1 <- ggplot() +
 geom_sf(data = firesheds.biv, 
         aes(fill = bi_class), 
         color = NA, linewidth = 0,
         show.legend = F) +
 bi_scale_fill(pal = "DkViolet", dim = 3) +
 theme_void() +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.15,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.25,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)

# create the legend
legend <- bi_legend(
 pal = "DkViolet",
 dim = 3,
 xlab = "Future Fire ",
 ylab = "Future Aspen ",
 pad_width = 1.2,
 size = 9)

# merge the plot and legend
p1.1 <- ggdraw() +
 draw_plot(p1) +
 draw_plot(legend, 0.17, .72, 0.25, 0.25)
p1.1

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_Bivar_FutureFireAspen.png')
ggsave(out_png, plot = p1.1, dpi = 500, width = 6, height = 8, bg = 'white')


###########################################
# Count occurrences of each bivariate class
bi_counts <- future.fire.bi %>%
 as_tibble() %>%
 count(bi_class)

# Bar chart of bivariate class distribution
counts_plot <- ggplot(bi_counts, aes(x = bi_class, y = n, fill = bi_class)) +
 geom_col(show.legend = FALSE) +
 bi_scale_fill(pal = "DkViolet", dim = 3) +  # Use the same bivariate color scheme
 theme_minimal() +
 labs(x = "Bivariate Class", y = "Cell Count")
counts_plot

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_counts.png')
ggsave(out_png, plot = counts_plot, height = 4, width = 4, bg = 'white')


#############
# WUI summary
wui_summary <- future.fire.bi %>%
 as_tibble() %>%
 group_by(bi_class) %>%
 summarize(
  WUI_Interface = mean(wui3, na.rm = TRUE),
  WUI_Intermix = mean(wui4, na.rm = TRUE),
  Distance_to_WUI = mean(wui_dist_mean, na.rm = TRUE)
 ) %>%
 pivot_longer(cols = c(WUI_Interface, WUI_Intermix), 
              names_to = "WUI_Class", values_to = "Percent_Cover")
glimpse(wui_summary)

# plot it
wui_barplot <- ggplot(wui_summary, aes(x = bi_class, y = Percent_Cover, fill = WUI_Class)) +
 geom_col(position = "dodge", color = "black") +  # Grouped bars with outlines
 scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +  # Distinct colors
 theme_minimal(base_size = 12) +
 labs(x = "Bivariate Class", y = "Mean % Cover", fill = "WUI Class") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))
wui_barplot

# distance to wui
wui_boxplot <- ggplot(future.fire.bi, aes(x = bi_class, y = wui_dist_mean, fill = bi_class)) +
 geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxes
 geom_jitter(alpha = 0.2, width = 0.2, size = 0.8) +  # Light scatter points
 bi_scale_fill(pal = "DkViolet", dim = 3) +  # Keep consistent with bivariate map
 theme_minimal(base_size = 12) +
 labs(x = "Bivariate Class", y = "Distance to WUI (m)") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))
wui_boxplot


# create a combined plot
# Arrange the two summary plots in a vertical stack (right side)
# Merge with the bivariate map (left side)
final_plot <- plot_grid(
 p1.1,        # Bivariate map
 counts_plot,  # Stacked plots
 ncol = 2,    # Two-column layout
 labels = c("A", "B"), 
 rel_widths = c(2, 1.3)
)
final_plot

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_wInsets.png')
ggsave(out_png, plot = final_plot, dpi = 500, bg = 'white')

