
library(tidyverse)
library(sf)

# Animation libraries
library(gganimate)
library(gifski)
library(transformr)

# Data

getwd()

frp <- st_read('data/spatial/raw/VIIRS/DL_FIRE_SV-C2_342931/fire_archive_SV-C2_342931.shp')
glimpse(frp)

# Create static map
static <- ggplot() +
  # geom_sf(data = frp, fill="grey90", size = 0.8, col = NA) +
  geom_sf(data = frp, aes(color=FRP, size=FRP), alpha=0.4) +
  scale_color_viridis_c(option = "plasma", trans = "log10") + 
                        # limits = c(1864,2019), 
                        # breaks = c(1864, 1900, 1950, 1990, 2000, 2019),) +
  scale_size_continuous(trans="log10", range = c(0.2,1.4), guide=NULL) +
  guides(color = guide_colourbar(position="top", barwidth = 0.8, barheight = 12.5, ticks=F,
                                 label.position = "right", title.position = "left",
                                 label.theme = element_text(angle = 0, size=12))) +
  ggspatial::annotation_scale(height=unit(1.2, "mm")) +
  labs(
       color="FRP") +
  theme_void() +
  theme(
        legend.title = element_text(angle = 90, size=14))
static

# Animate it
animap <- static +
  transition_time(ACQ_DATE) +
  labs(caption = "Date: {frame_time}\nFrame {frame} of {nframes}") +
  gganimate::enter_recolor(fill = "#f0f5f9") +
  gganimate::shadow_mark(past = TRUE, alpha = 0.8) +
  theme(plot.caption = element_text(size=12))
# Grab number of frames
ndays <- length(unique((frp$ACQ_DATE))) + 25
gganimate::animate(
  animap, nframes = ndays, fps=8, 
  start_pause = 3, end_pause=10,
  renderer = gifski_renderer(),
  width=650,height=650)
# Save GIF
anim_save("frp_animation_mullen_fire.gif")

