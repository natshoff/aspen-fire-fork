# # Read in the forest type summaries (TreeMap 2016)
# 
# treemap <- read_csv("data/tabular/mod/TreeMap/frp_plot_foresttype.csv") %>%
#  mutate(
#   # Binary aspen presence
#   presence = if_else(ForTypName == "Aspen", 1, 0)
#  ) %>%
#  # Pivot to wider, with land cover as column names
#  pivot_wider(
#   names_from = ForTypName,
#   values_from = FTypePct,
#   values_fill = 0) %>%
#  select(FID,`Aspen`,FTypeAcres,presence) %>%
#  rename(fid = FID,
#         aspen_pct = `Aspen`,
#         aspen_acres = FTypeAcres) %>%
#  group_by(fid) %>%
#  summarize(aspen_pct = sum(aspen_pct),
#            aspen_pres = max(presence))
# head(treemap)
# 
# # Join to the FRP data frame
# frp <- frp %>% left_join(treemap,by="fid") %>%
#  mutate(aspen_pct = if_else(is.na(aspen_pct), 0, aspen_pct),
#         aspen_pres = if_else(is.na(aspen_pct), 0, 1))
# 
# rm(treemap)
# 
# # Read in the GEE reductions
# 
# # Fire weather (daily)
# gmet <- read_csv("data/tabular/mod/gee-exports/vnp14img_west_weather-daily.csv") %>%
#  select(fid,bi_max,vpd_max,vs_max) %>%
#  drop_na()
# # Fire weather (long-term)
# gmet_lt <- read_csv("data/tabular/mod/gee-exports/vnp14img_west_weather-lt.csv") %>%
#  select(fid,bi_p90,vpd_p90) %>%
#  drop_na()
# # Topography
# topo <- read_csv('data/tabular/mod/gee-exports/vnp14img_west_topo.csv') %>%
#  select(fid,elev,slope,roughness) %>%
#  drop_na()
# # Previous season EVI
# evi <- read_csv("data/tabular/mod/gee-exports/vnp14img_west_veg-evi.csv") %>%
#  select(fid,evi) %>%
#  drop_na()
# # LCMS modal land cover
# lc <- read_csv("data/tabular/mod/gee-exports/vnp14img_west_veg-lc.csv") %>%
#  select(fid,lcms) %>%
#  drop_na()
# # MODIS vegetated continuous fields (VCF)
# vcf <- read_csv("data/tabular/mod/gee-exports/vnp14img_west_veg-vcf.csv") %>%
#  select(fid, pct_notree_veg, pct_noveg, pct_tree) %>%
#  drop_na()
# 
# 
# # Bind back to the FRP observations
# frp <- frp %>%
#  left_join(gmet,by="fid") %>%
#  left_join(topo,by="fid") %>%
#  left_join(evi,by="fid") %>%
#  left_join(lc,by="fid") %>%
#  left_join(vcf,by="fid") %>%
#  select(
#   fid,ACQ_DATE,DAYNIGHT,BRIGHTNESS,FRP,CONFIDENCE,
#   Event_ID,Ig_Date,BurnBndAc,
#   mtbs_aspen_ac,mtbs_aspen_pct,aspen_pct,
#   bi_max,vpd_max,vs_max,
#   cbh_mean,cbh_stdev,cbd_mean,cbd_stdev,
#   evi,lcms,pct_notree_veg,pct_noveg,pct_tree,
#   elev,slope,roughness) %>%
#  rename(
#   id = fid,
#   mtbs_id = Event_ID,
#   mtbs_date = Ig_Date,
#   acq_date = ACQ_DATE,
#   daynight = DAYNIGHT,
#   brightness = BRIGHTNESS,
#   confidence = CONFIDENCE,
#   frp = FRP,
#   cbh_mn = cbh_mean,
#   cbh_sd = cbh_stdev,
#   cbd_mn = cbd_mean,
#   cbd_sd = cbd_stdev,
#   evi_mn = evi,
#   lcms_mode = lcms
#  ) %>%
#  mutate(mtbs_date = as.Date(mtbs_date)) %>%
#  st_transform(st_crs(mtbs)) %>%
#  # Keep only complete cases for now ... 
#  # Likely need to investigate why there are so many NAs ...
#  drop_na()
# glimpse(frp)
# 
# # Check on completeness
# 
# 
# # Tidy up!
# rm(mtbs,gmet,gmet_lt,topo,evi,lc,vcf)
# gc()
# 
# # Write out the new file
# st_write(frp,"data/spatial/mod/vnp14img_west_spatial_w_attr.gpkg",delete_dsn=TRUE,append=FALSE)

