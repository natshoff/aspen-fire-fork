
# Load the libraries
library(tidyverse)
library(sf)
library(biscale)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/"

# Load the data (future fire, future aspen)
future.fire <- st_read(paste0(projdir, 'Aim3/data/spatial/mod/future_fire_grid.gpkg'))
glimpse(future.fire)

# create the bivariate classes
future.fire.bi <- bi_class(future.fire, x = En_NFire, y = f_aspen_p90, style = "quantile", dim = 3)

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

