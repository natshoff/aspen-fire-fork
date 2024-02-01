
library(tidyverse)
library(lubridate)
library(sf)

# Load the grids, keep the elevation attribute
grid <- st_read('data/spatial/mod/boundaries/spatial_block_grid_50km2.gpkg') %>%
  select(grid_id,elevation_mn,treemap_sum)
glimpse(grid)

# Load the phenology by spatial block grid
phenology <- read_csv('data/tabular/mod/phenology/viirs_phenology_by_grid_in_aspen.csv') %>%
  rename(fid = `system:index`) %>%
  mutate(grid_id = str_sub(fid, -7),
         year = as.integer(str_sub(fid, 1, 4))) %>%
  select(-c(fid,label,count,.geo)) %>%
  # Create DOY versions of the phenology metrics for plotting and statistical analysis
  mutate(
    season_length = Growing_Season_Length_1_median,
    doy_midgreenup = yday(Date_Mid_Greenup_Phase_1_median),
    doy_midsenescence = yday(Date_Mid_Senescence_Phase_1_median),
    doy_green_dec = yday(Onset_Greenness_Decrease_1_median),
    doy_green_inc = yday(Onset_Greenness_Increase_1_median),
    doy_green_max = yday(Onset_Greenness_Maximum_1_median),
    doy_green_min = yday(Onset_Greenness_Minimum_1_median)
  ) %>%
  left_join(grid%>%as_tibble(), by="grid_id") %>%
  select(grid_id,year,elevation_mn,season_length,doy_midgreenup,doy_midsenescence,
         doy_green_dec,doy_green_inc,doy_green_max,doy_green_min)
glimpse(phenology)

# Reshape the data frame for plotting
# Add a statistics column for legend (could do this for F1 too ...)
phenology.m <- reshape2::melt(
  phenology, id.vars = c("grid_id","year","elevation_mn"), 
  measure.vars = c("season_length","doy_midgreenup","doy_midsenescence",
                   "doy_green_dec","doy_green_inc","doy_green_max","doy_green_min"), 
  variable.name = "metric",
  value.name = "doy"
) 
head(phenology.m)

# Create a growing season length plot by grid ID
ggplot(phenology, aes(x=factor(year), y=season_length)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Growing Season Length",
       x = "Year",
       y = "Growing Season Length",
       fill = "Year") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Growing season length as a function of elevation
ggplot(data=phenology, aes(x=elevation_mn,y=doy_green_max)) +
    geom_smooth(method="lm",colour="gray40", fill="gray70", size=0.8) +
    geom_point(size=0.4, color="gray10", fill=NA) +
    facet_wrap(~factor(year), scales="free_x") +
    theme_light(12) +
    labs(title="Onset of Maximum Greenness by Elevation",
         x="Elevation",y="Day-of-Year")

# DOY plots

# Facet wrap plot

# Reshape from wide to long format
phenology.l <- phenology %>%
  pivot_longer(
    cols = c(doy_midgreenup,doy_midsenescence,doy_green_dec,
             doy_green_inc,doy_green_max,doy_green_min),
    names_to = "metric",
    values_to = "doy"
  )
head(phenology.l)
glimpse(phenology.l)

ordered_metrics <- c(
  "doy_green_inc","doy_midgreenup","doy_green_max",
  "doy_green_dec","doy_midsenescence","doy_green_min") 

metric_names <- c(
  doy_green_inc = "Onset Greenness Increase",
  doy_midgreenup = "Mid Greenup Phase",
  doy_green_max = "Onset Greenness Maximum",
  doy_green_dec = "Onset Greenness Decrease",
  doy_midsenescence = "Mid Senescence Phase",
  doy_green_min = "Onset Greenness Minimum"
)

phenology.l$metric <- factor(
  phenology.l$metric, 
  levels = ordered_metrics, 
  labels = metric_names
)

# Create the plot with facet_wrap
ggplot(phenology.l, aes(x=factor(year), y=doy)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "fixed") +
  theme_minimal() +
  labs(title = "DOY Phenology Metrics (2013-2022)",
       x = "Year",
       y = "Day of Year",
       fill = "Year") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Stacked box plots
ggplot(phenology.l, aes(x=factor(year), y=doy, fill=metric)) +
  geom_boxplot(width=1) +
  theme_minimal() +
  labs(x = "Year",
       y = "Day of Year",
       fill = "Metric") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# Relationship with elevation
ggplot(data=phenology.l, aes(x=elevation_mn,y=doy,color=metric)) +
  geom_point(size=1.2) +
  theme_light(12)

# Calculate the average metrics to help in determining the optimal windows
phenology_by_year <- phenology %>%
  select(year, starts_with('doy_')) %>%
  group_by(year) %>%
  summarize_all(., ~ as.integer(mean(., na.rm = TRUE)))
glimpse(phenology_by_year)

print("Metric averages across blocks and years: ")
summary(phenology_by_year)

# Grab the 2019 summary across all blocks
summary(phenology%>%filter(year==2019))

