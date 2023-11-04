
library(tidyverse)
library(sf)

getwd()

proj <- st_crs("ESRI:102039")

# Tidy the model data frame (daily FRP observations)

# Spatial data
frp.sp <- st_read("data/spatial/mod/vnp14img_west_spatial_w_attr.gpkg") %>%
 # tidy the frame, select model coefficients
 select(
  id,mtbs_id,acq_date,frp,aspen_pct,cbd_mn,cbh_mn,
  vs_max,bi_max,vpd_max,elev,slope,roughness,
  evi_mn,pct_tree,pct_notree_veg,pct_noveg,
  geom
 ) %>%
 # Convert back to centroid (points)
 st_centroid() %>%
 st_transform(proj)

# Tabular format
frp.df <- frp.sp %>% st_set_geometry(NULL)
# Center and scale the model coefficients
frp.df <- frp.df %>%
 mutate_at(c("cbd_mn","cbh_mn","evi_mn",
             "vs_max","bi_max","vpd_max",
             "elev","slope","roughness",
             "pct_tree","pct_notree_veg","pct_noveg"), 
           ~(scale(.) %>% as.vector))
glimpse(frp.df)


# Calculate the IPW for aspen percent cover
# Do this for each fire event individually

dfs <- list()
fires <- unique(frp.df$mtbs_id)
for (i in 1:length(fires)) {
 print(fires[i])
 
 # Get the fire data frame, y and x variables
 df <- frp.df %>% filter(mtbs_id == fires[i]) %>% drop_na()
 print(dim(df))
 
 num <- 0
 den <- 0
 
 # Calculate the IPW for the fire
 # Calculate the numerator (expected distribution of aspen %)
 model_num <- lm(aspen_pct ~ 1, data = df)
 num <- dnorm(df$aspen_pct, predict(model_num), sd(model_num$residuals))
 
 # Now calculate the denominator (exp. distribution regressed against confounders)
 model_den <- lm(
  aspen_pct ~ cbd_mn + cbh_mn + vs_max + bi_max + vpd_max + 
   elev + slope + roughness + evi_mn + 
   pct_tree + pct_notree_veg + pct_noveg, data=df)
 den <- dnorm(df$aspen_pct, predict(model_den), sd(model_den$residuals))
 
 # Add the IPW column to the data frame
 df <- df %>% mutate(ipw = num / den)
 
 dfs[[i]] <- df
}


