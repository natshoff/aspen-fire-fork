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
library(ggcorrplot)
library(RColorBrewer)
library(viridis)

# Install and load the pracma package if you haven't
if(!require(pracma)) {
  install.packages("pracma")
  library(pracma)
}

setwd('~/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1')

# boundary data
srme <- st_read("data/spatial/mod/boundaries/us_eco_l3_srme.gpkg")
wrnf <- st_read("data/case-study/wrnf_boundary_100m.gpkg")
blocks <- st_read("data/spatial/mod/boundaries/spatial_block_grid_w_attr.gpkg") %>%
  mutate(Dormancy_1 = as.Date(Dormancy_1, origin = "1970-01-01"),
         Greenup_1 = as.Date(Greenup_1, origin = "1970-01-01"),
         Maturity_1 = as.Date(Maturity_1, origin = "1970-01-01"),
         MidGreendown_1 = as.Date(MidGreendown_1, origin = "1970-01-01"),
         MidGreenup_1 = as.Date(MidGreenup_1, origin = "1970-01-01"),
         Peak_1 = as.Date(Peak_1, origin = "1970-01-01"),
         Senescence_1 = as.Date(Senescence_1, origin = "1970-01-01")) %>%
  rename(elevation_mn = mean)

# Spectral Response
ts <- read_csv('data/tabular/mod/results/s2l2a_w_indices-time_series-y19.csv')

# Accuracy Metrics
# source('code/R/accmeas.R')
# source('code/R/accmeas_sensitivity.R')
accmeas <- read_csv('data/tabular/mod/results/accmeas_prop.csv')
accmeas.s <- read_csv('data/tabular/mod/results/sensitivity/accmeas_f1best_sensitivity.csv')

# Feature Importance
ftr_imp <- read_csv('data/tabular/mod/results/srme_skcv_ftr_imps.csv')

# Landscape Patch Analysis

# Test data (Sentinel-based aspen map)

# 10-meter
ls_ptch_10m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_patch_metrics.csv") %>%
  mutate(source = "Aspen10m")
ls_cls_10m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_class_metrics.csv") %>%
  mutate(source = "Aspen10m")

# 30-meter maximum
ls_ptch_30m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_max30m_patch_metrics.csv") %>%
  mutate(source = "Aspen30m")
ls_cls_30m <- read.csv("data/tabular/mod/results/ls_metrics/aspen_prob_10m_binOpt_max30m_class_metrics.csv") %>%
  mutate(source = "Aspen30m")

# Reference data (LFEVT, TreeMap, ITSP)

# LANDFIRE
ls_ptch_evt <- read.csv("data/tabular/mod/results/ls_metrics/lc16_evt_srme_aspen_r01_utm_patch_metrics.csv") %>%
  mutate(source = "LFEVT")
ls_cls_evt <- read.csv("data/tabular/mod/results/ls_metrics/lc16_evt_srme_aspen_r01_utm_class_metrics.csv") %>%
  mutate(source = "LFEVT")

# ITSP (basal area > 0)
ls_ptch_itsp <- read.csv("data/tabular/mod/results/ls_metrics/itsp_aspen_srme_r1__patch_metrics.csv") %>%
  mutate(source = "ITSP")
ls_cls_itsp <- read.csv("data/tabular/mod/results/ls_metrics/itsp_aspen_srme_r1__class_metrics.csv") %>%
  mutate(source = "ITSP")

# USFS TreeMap
ls_ptch_treemap <- read.csv("data/tabular/mod/results/ls_metrics/treemap16_spcd_746_int_patch_metrics.csv") %>%
  mutate(source = "TreeMap")
ls_cls_treemap <- read.csv("data/tabular/mod/results/ls_metrics/treemap16_spcd_746_int_class_metrics.csv") %>%
  mutate(source = "TreeMap")

# Create one data frame for the landscape metrics

# Patch metrics
patch_metrics <- ls_ptch_evt %>%
  bind_rows(ls_ptch_itsp) %>%
  bind_rows(ls_ptch_treemap) %>%
  bind_rows(ls_ptch_10m) %>%
  bind_rows(ls_ptch_30m)

# Class metrics
class_metrics <- ls_cls_evt %>%
  bind_rows(ls_cls_itsp) %>%
  bind_rows(ls_cls_treemap) %>%
  bind_rows(ls_cls_10m) %>%
  bind_rows(ls_cls_30m)

# Clean up !
rm(
  ls_ptch_evt,ls_cls_evt,ls_ptch_itsp,ls_cls_itsp,
  ls_ptch_treemap,ls_cls_treemap,
  ls_ptch_10m,ls_ptch_30m,ls_cls_10m,ls_cls_30m
)


# Agreement assessment

# Global

ref.srme <- read_csv("data/tabular/mod/results/global_accmeas_multi_blocks_full_srme.csv") %>%
  mutate(region = "SRME")
ref.wrnf <- read_csv("data/tabular/mod/results/global_accmeas_multi_blocks_full_wrnf.csv") %>%
  mutate(region = "WRNF")
ref.global <- bind_rows(ref.srme,ref.wrnf) %>%
  mutate(
    source = if_else(source=="lc16_evt_srme_aspen_r01_utm_match_wrnf", "LFEVT", source),
    source = if_else(source=="treemap16_spcd_746_int_match_wrnf", "TreeMap", source),
    source = if_else(source=="itsp_aspen_srme_r1__match_wrnf", "ITSP", source),
    source = if_else(source=="lc16_evt_srme_aspen_r01_utm_match_srme", "LFEVT", source),
    source = if_else(source=="treemap16_spcd_746_int_match_srme", "TreeMap", source),
    source = if_else(source=="itsp_aspen_srme_r1__match_srme", "ITSP", source),
    # Calculate the F1 Score
    f1 = 2 * (prec * rec) / (prec + rec)
  )
glimpse(ref.global)
rm(ref.srme,ref.wrnf)

# Focal

lf.focal.srme <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_lc16_evt_srme_aspen_r01_utm_match_srme.csv') %>%
  mutate(region = "SRME")
lf.focal.wrnf <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_lc16_evt_srme_aspen_r01_utm_match_wrnf.csv') %>%
  mutate(region = "WRNF")

itsp.focal.srme <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_itsp_aspen_srme_r1__match_srme.csv') %>%
  mutate(region = "SRME")
itsp.focal.wrnf <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_itsp_aspen_srme_r1__match_wrnf.csv') %>%
  mutate(region = "WRNF")

treemap.focal.srme <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_treemap16_spcd_746_int_match_wrnf.csv') %>%
  mutate(region = "SRME")
treemap.focal.wrnf <- read_csv('data/tabular/mod/results/focal_accmeas_multi_blocks_treemap16_spcd_746_int_match_wrnf.csv') %>%
  mutate(region = "WRNF")

# Bind them together
ref.focal <- bind_rows(lf.focal.srme,lf.focal.wrnf,itsp.focal.srme,itsp.focal.wrnf,treemap.focal.srme,treemap.focal.wrnf)
glimpse(ref.focal)

# Clean up
rm(lf.focal.srme,lf.focal.wrnf,itsp.focal.srme,itsp.focal.wrnf,treemap.focal.srme,treemap.focal.wrnf)
