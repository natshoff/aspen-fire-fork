
# Load the libraries
library(tidyverse)
library(sf)
library(biscale)
library(cowplot)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/"

# Load the data (future fire, future aspen)
future.fire <- st_read(paste0(projdir, 'Aim3/data/spatial/mod/future_fire_grid_trend_aspen.gpkg'))
glimpse(future.fire)

#######################
# aspen suitability map
p1 <- ggplot(future.fire) +
 geom_sf(aes(fill = f_aspen_mn), color = NA) +  # No border for smooth visualization
 scale_fill_viridis(option = "viridis", na.value = "gray80", name = "Aspen Suitability") +
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
# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_FutureAspen_Grid.png')
ggsave(out_png, plot = p1, dpi = 500, bg = 'white')


############################################
# Bivariate map: Future Fire -> Future Aspen
# create the bivariate classes
future.fire.bi <- bi_class(future.fire, x = trend_count_sc, y = f_aspen_mn_sc, style = "quantile", dim = 3)

# create the bivariate chloropleth map
p1 <- ggplot() +
 geom_sf(data = future.fire.bi, 
         aes(fill = bi_class), 
         color = NA, linewidth = 0,
         show.legend = F) +
 bi_scale_fill(pal = "DkViolet", dim = 3) +
 theme_void()

# create the legend
legend <- bi_legend(
 pal = "DkViolet",
 dim = 3,
 xlab = "Future Fire ",
 ylab = "Future Aspen ",
 pad_width = 1.5,
 size = 12)

# merge the plot and legend
p1.1 <- ggdraw() +
 draw_plot(p1) +
 draw_plot(legend, 0.17, .72, 0.25, 0.25)
p1.1

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen.png')
ggsave(out_png, plot = p1.1, dpi = 500, bg = 'white')


# Count occurrences of each bivariate class
bi_counts <- future.fire.bi %>%
 as_tibble() %>%
 count(bi_class)

# Bar chart of bivariate class distribution
counts_plot <- ggplot(bi_counts, aes(x = bi_class, y = n, fill = bi_class)) +
 geom_col(show.legend = FALSE) +
 bi_scale_fill(pal = "DkViolet", dim = 3) +  # Use the same bivariate color scheme
 theme_bw() +
 labs(x = "Bivariate Class", y = "Cell Count")
counts_plot

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_counts.png')
ggsave(out_png, plot = counts_plot, height = 6, width = 4, bg = 'white')

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

