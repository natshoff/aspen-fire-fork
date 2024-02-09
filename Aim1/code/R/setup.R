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

#############
# Functions #
#############

# Function to calculate accuracy and F1 score
calculate_metrics <- function(df, true_label, pred_label) {
  tp <- sum((df[[true_label]] == 1) & (df[[pred_label]] == 1))
  tn <- sum((df[[true_label]] == 0) & (df[[pred_label]] == 0))
  fp <- sum((df[[true_label]] == 0) & (df[[pred_label]] == 1))
  fn <- sum((df[[true_label]] == 1) & (df[[pred_label]] == 0))
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(list(accuracy = accuracy, f1_score = f1_score))
}


########
# Data #
########

# Spatial data
srme <- st_read("data/spatial/raw/boundaries/us_eco_l3_srme.gpkg")
wrnf <- st_read("data/spatial/raw/boundaries/wrnf_boundary.gpkg")
blocks <- st_read("data/spatial/mod/boundaries/spatial_block_grid_50km2_count_s2.gpkg")


############################
# Landscape Patch Analysis #
############################

# Southern Rockies

ls_dir <- 'data/tabular/mod/results/ls_metrics/patch'
ls_patch_files.sr <- list.files(
  path = ls_dir, 
  pattern = "srme", 
  full.names = TRUE)

ls_dir <- 'data/tabular/mod/results/ls_metrics/class'
# Southern Rockies
ls_class_files.sr <- list.files(
  path = ls_dir, 
  pattern = "srme", 
  full.names = TRUE)

# Load into a data frame
patch.sr <- bind_rows(map(ls_patch_files.sr, read_csv))
class.sr <- bind_rows(map(ls_class_files.sr, read_csv))


# White River NF

ls_dir <- 'data/tabular/mod/results/ls_metrics/patch'
# Southern Rockies
ls_patch_files.wr <- list.files(
  path = ls_dir, 
  pattern = "wrnf", 
  full.names = TRUE)

ls_dir <- 'data/tabular/mod/results/ls_metrics/class'
# Southern Rockies
ls_class_files.wr <- list.files(
  path = ls_dir, 
  pattern = "wrnf", 
  full.names = TRUE)

# Load into a data frame
patch.wr <- bind_rows(map(ls_patch_files.wr, read_csv))
class.wr <- bind_rows(map(ls_class_files.wr, read_csv))


########################
# Agreement assessment #
########################

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

