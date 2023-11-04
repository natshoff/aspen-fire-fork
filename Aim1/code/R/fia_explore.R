
###################

devtools::install_github('hunter-stanke/rFIA')

require(rFIA)
require(tidyverse)
require(sf)
library(DBI)

###################

setwd("/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1/data/tabular/raw/FIA")

##############################################################################
## Title: How to save tables from a SQLite database as individual CSVs`
## Author: Karin Kralicek (karin.kralicek@usda.gov)
## Date: 04/11/2022
##
## About: 
## - Work-around for issue with Datamart downloads for CSVs, which are giving a
##   404 error at present.
###############################################################################

unlistFIA <- function(st){
  # Set path to where the SQLite is located
  path_db <- paste0(getwd(), "/SQLite_FIADB_", st, "/")
  
  # Set path to where you want the CSVs to appear
  path_out = paste0(path_db, "tables/")
  create(path_out)
  
  # Take a reference to the db (again, using CO as the example)
  con <- dbConnect(RSQLite::SQLite(), 
                   paste0(path_db, "FIADB_CO.db"))
  
  # grab all the table names
  # - alternatively, subset this to only those tables you want
  #   e.g. `db_table_names <- c("SURVEY", "PLOT")`
  db_table_names <- dbListTables(con)
  
  # iterate through the table names and write out a csv for each table in the db
  lapply(db_table_names, function(x) {
    write.csv(dbReadTable(conn = con, x),
              file = paste0(path_db, "CO_", x, ".csv"))
  })
  
  # close connection
  dbDisconnect(con)
}

# Run it ...

#################################################################################

# Now work with the FIA data tables
st <- "CO"
f <- paste(getwd(), "/SQLite_FIADB_", st, "/tables/", sep="")
# Read the FIA data
db <- rFIA::readFIA(f)
# Read the spatial bounds (population)
hucs <- st_read("../../../spatial/raw/boundaries/srm_huc10_blue.shp")

# Clip to our spatial domain
dbc <- clipFIA(db, mask = hucs, mostRecent = FALSE)
dbc_spat <- tpa(dbc, polys = hucs, bySpecies = TRUE, returnSpatial = TRUE)
dbc.potr <- dbc_spat %>% filter(SCIENTIFIC_NAME == "Populus tremuloides")
# produce a sample map
ggplot() +
  geom_sf(data=dbc.potr, aes(fill=BAA)) +
  scale_fill_viridis_c(option="inferno", trans="log10")





