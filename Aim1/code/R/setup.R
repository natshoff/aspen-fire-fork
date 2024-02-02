############
# LIBRARIES
############

library(tidyverse)
library(sf)
library(lubridate)
library(scales)
library(ggfortify)
library(gridExtra)
library(ggthemes)
library(ggpubr)
library(randomForest)
library(rfUtilities)
library(caret)
library(ROCR)
library(RColorBrewer)
library(viridis)

# Install and load the pracma package if you haven't
if(!require(pracma)) {
  install.packages("pracma")
  library(pracma)
}

setwd('~/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1')

# Spatial data
srme <- st_read("data/spatial/raw/boundaries/us_eco_l3_srme.gpkg")
wrnf <- st_read("data/spatial/raw/boundaries/wrnf_boundary.gpkg")
blocks <- st_read("data/spatial/mod/boundaries/spatial_block_grid_50km2.gpkg")

# Accuracy Metrics
accmeas <- read_csv('data/tabular/mod/results/best_model/southern_rockies_accmeas.csv',
                    show_col_types = FALSE)

# Feature Importance
ftr_imp <- read_csv('data/tabular/mod/results/best_model/southern_rockies_feature_imps.csv',
                    show_col_types = FALSE)


# # Landscape Patch Analysis
# 
# # Test data (Sentinel-based aspen map)
# 
# # 10-meter
# ls_ptch_10m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_patch_metrics.csv") %>%
#   mutate(source = "Aspen10m")
# ls_cls_10m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_class_metrics.csv") %>%
#   mutate(source = "Aspen10m")
# 
# # 30-meter maximum
# ls_ptch_30m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_max30m_patch_metrics.csv") %>%
#   mutate(source = "Aspen30m")
# ls_cls_30m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_max30m_class_metrics.csv") %>%
#   mutate(source = "Aspen30m")
# 
# # Reference data (LFEVT, TreeMap, ITSP)
# 
# # LANDFIRE
# ls_ptch_evt <- read.csv("data/tabular/mod/results/ls_metrics/lc16_evt_srme_aspen_r01_utm_patch_metrics.csv") %>%
#   mutate(source = "LFEVT")
# ls_cls_evt <- read.csv("data/tabular/mod/results/ls_metrics/lc16_evt_srme_aspen_r01_utm_class_metrics.csv") %>%
#   mutate(source = "LFEVT")
# 
# # ITSP (basal area > 0)
# ls_ptch_itsp <- read.csv("data/tabular/mod/results/ls_metrics/itsp_aspen_srme_r1__patch_metrics.csv") %>%
#   mutate(source = "ITSP")
# ls_cls_itsp <- read.csv("data/tabular/mod/results/ls_metrics/itsp_aspen_srme_r1__class_metrics.csv") %>%
#   mutate(source = "ITSP")
# 
# # USFS TreeMap
# ls_ptch_treemap <- read.csv("data/tabular/mod/results/ls_metrics/treemap16_spcd_746_int_patch_metrics.csv") %>%
#   mutate(source = "TreeMap")
# ls_cls_treemap <- read.csv("data/tabular/mod/results/ls_metrics/treemap16_spcd_746_int_class_metrics.csv") %>%
#   mutate(source = "TreeMap")
# 
# # Create one data frame for the landscape metrics
# 
# # Patch metrics
# patch_metrics <- ls_ptch_evt %>%
#   bind_rows(ls_ptch_itsp) %>%
#   bind_rows(ls_ptch_treemap) %>%
#   bind_rows(ls_ptch_10m) %>%
#   bind_rows(ls_ptch_30m)
# 
# # Class metrics
# class_metrics <- ls_cls_evt %>%
#   bind_rows(ls_cls_itsp) %>%
#   bind_rows(ls_cls_treemap) %>%
#   bind_rows(ls_cls_10m) %>%
#   bind_rows(ls_cls_30m)
# 
# # Clean up !
# rm(
#   ls_ptch_evt,ls_cls_evt,ls_ptch_itsp,ls_cls_itsp,
#   ls_ptch_treemap,ls_cls_treemap,
#   ls_ptch_10m,ls_ptch_30m,ls_cls_10m,ls_cls_30m
# )
# 
# 
# Agreement assessment

# Global

aggr_dir <- 'data/tabular/mod/results/agreement'
file_list.srme <- list.files(path = aggr_dir, pattern = "\\srme_10m.csv$", full.names = TRUE)
file_list.wrnf <- list.files(path = aggr_dir, pattern = "\\wrnf_10m.csv$", full.names = TRUE)

# Read and merge the CSV files for each region

aggr.sr <- read_csv('data/tabular/mod/results/agreement/global_accmeas_multi_blocks_full_srme.csv') %>%
  rename(f1 = fi) %>%
  mutate(
    region = "Southern Rockies",
    source = if_else(source=="lc16_evt_200_bin_srme_10m", "LANDFIRE EVT", source),
    source = if_else(source=="usfs_itsp_aspen_ba_gt10_srme_10m", "USFS ITSP", source),
    source = if_else(source=="usfs_treemap16_bin_srme_10m", "USFS TreeMap", source))
glimpse(aggr.sr)

# White River NF
aggr.wr <- read_csv('data/tabular/mod/results/agreement/global_accmeas_multi_blocks_full_wrnf.csv') %>%
  mutate(
    region = "White River NF",
    source = if_else(source=="lc16_evt_200_bin_wrnf_10m", "LANDFIRE EVT", source),
    source = if_else(source=="usfs_itsp_aspen_ba_gt10_wrnf_10m", "USFS ITSP", source),
    source = if_else(source=="usfs_treemap16_bin_wrnf_10m", "USFS TreeMap", source),
    # Calculate the F1 Score
    f1 = 2 * (prec * rec) / (prec + rec)) %>%
  rename(precision = prec,
         recall = rec)
glimpse(aggr.wr)

