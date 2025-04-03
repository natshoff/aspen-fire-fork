
# shiny app for SOM analysis

# libraries
library(shiny)
library(kohonen)
library(tidyverse)
library(sf)
library(scales)
library(ggradar)

# # set the working directory for the app
# projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"
# setwd(paste0(projdir,"code/R/app/"))

# Load firesheds
firesheds <- st_read('data/srm_firesheds_model_data.gpkg')

# Labels for the input variables
var_labels <- c(
 trend_count = "Future fire occurrence (trend)",
 delta585 = "Aspen suitability change (SSP585)",
 aspen_forest = "Aspen forest proportion",
 aspen10_ha = "Aspen absolute area (ha)",
 n_patches = "Aspen number of patches",
 patch_den = "Aspen patch density",
 big_patch = "Largest patch size",
 combust_sum = "Combustible mass",
 pop_den = "Population density",
 pop_n = "Population count",
 wui = "WUI %",
 wui_dist = "Distance to WUI",
 burned_pct_c = "Burned % Since 1984 (cumulative)",
 whp_p90 = "Wildfire Hazard Potential (90th)",
 hui_p90 = "Housing Unit Impact (90th)",
 cfl_p90 = "Conditional Flame Length (90th)",
 forest_pct = "Proportion forested area",
 lf_canopy = "Canopy cover (%)",
 aspen = "Aspen (LANDFIRE)",
 douglas_fir = "Douglas-fir (LANDFIRE)",
 lodgepole = "Lodgepole pine (LANDFIRE)",
 gambel_oak = "Gambel oak (LANDFIRE)",
 sagebrush = "Sagebrush (LANDFIRE)",
 pinon_juniper = "Pi\u00f1on-juniper (LANDFIRE)",
 ponderosa = "Ponderosa pine (LANDFIRE)",
 spruce_fir = "Spruce-fir (LANDFIRE)",
 white_fir = "White fir (LANDFIRE)"
)

# Prepare input data
X <- firesheds %>%
 st_drop_geometry() %>%
 rename(
  n_patches = number_of_patches, 
  patch_den = patch_density, 
  big_patch = largest_patch_index,
  pop_den = pop_density_max, 
  pop_n = pop_count_sum, 
  wui_dist = wui_dist_mean,
  lf_canopy = forest_cc_mean
 ) %>%
 mutate(aspen10_ha = aspen10_pixn * 0.01,
        forest_ha = forest_pixels * 0.09,
        aspen_forest = aspen10_ha / forest_ha,
        wui = wui1 + wui2 + wui3 + wui4) %>%
 select(names(var_labels)) %>%
 mutate(across(everything(), ~ replace_na(.x, 0))) %>%
 mutate(across(everything(), ~ as.numeric(scale(., center = TRUE))))

cluster_colors <- RColorBrewer::brewer.pal(6, "Accent")
names(cluster_colors) <- paste0("Cluster ", 1:6)

# leaflet setup
# UI
ui <- fluidPage(
 titlePanel("Aspen Firesheds-SOM Explorer"),
 
 fluidRow(
  column(
   width = 3,
   h4("Set Variable Weights"),
   actionButton("run_som", "Update SOM"),
   br(), br(),
   div(style = "overflow-y: auto; max-height: 550px; padding-right: 10px;",
       uiOutput("weight_sliders")
   ),
   numericInput("xdim", "Grid X dimension", 20, 5, 50),
   numericInput("ydim", "Grid Y dimension", 20, 5, 50)
  ),
  
  column(
   width = 9,
   fluidRow(
    column(6, plotOutput("som_map", height = "450px")),
    column(6, plotOutput("radar_plot", height = "450px"))
   )
  )
 )
)

server <- function(input, output, session) {
 vars <- colnames(X)
 
 output$weight_sliders <- renderUI({
  tagList(
   lapply(vars, function(v) {
    sliderInput(paste0("weight_", v), label = var_labels[v], min = 0, max = 5, value = 1, step = 0.1)
   })
  )
 })
 
 weights <- reactive({
  setNames(sapply(vars, function(v) input[[paste0("weight_", v)]]), vars)
 })
 
 som_results <- eventReactive(input$run_som, {
  X_weighted <- X * weights()
  som_grid <- somgrid(xdim = input$xdim, ydim = input$ydim, topo = "hexagonal")
  som_model <- som(as.matrix(X_weighted), grid = som_grid, rlen = 1000, alpha = c(0.1, 0.02))
  
  clusters <- cutree(hclust(dist(som_model$codes[[1]]), method = "ward.D2"), k = 6)
  cluster_assignments <- clusters[som_model$unit.classif]
  
  updated_firesheds <- firesheds
  updated_firesheds$som_cluster <- factor(paste0("Cluster ", cluster_assignments))
  
  # Cluster means for radar plot
  X_clustered <- X_weighted %>% mutate(cluster = updated_firesheds$som_cluster)
  cl_means <- X_clustered %>%
   group_by(cluster) %>%
   summarise(across(everything(), mean)) %>%
   mutate(across(-cluster, rescale)) %>%
   rename(group = cluster)
  
  list(map = updated_firesheds, radar = cl_means)
 })
 
 output$som_map <- renderPlot({
  ggplot(som_results()$map) +
   geom_sf(aes(fill = som_cluster), color = NA) +
   scale_fill_manual(values = cluster_colors, name = "SOM Cluster") +
   theme_void() +
   theme(legend.position = "bottom")
 })
 
 output$radar_plot <- renderPlot({
  req(som_results())
  ggradar(som_results()$radar,
          grid.min = 0, grid.mid = 0.5, grid.max = 1,
          group.line.width = 0.8, group.point.size = 1.5,
          axis.label.size = 3.5, grid.label.size = 4,
          group.colours = cluster_colors,
          fill = TRUE, fill.alpha = 0.2,
          legend.position = "none",
          font.radar = "sans")
 })
}

shinyApp(ui, server)

