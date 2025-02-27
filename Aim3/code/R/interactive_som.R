library(dplyr)
library(DT)

server <- function(input, output, session) {
 
 # Reactive weighted dataset
 weighted_data <- reactive({
  variable_weights <- c(
   wui1 = input$wui_weight, wui2 = input$wui_weight, 
   wui3 = input$wui_weight, wui4 = input$wui_weight,
   pop_density_mean = input$pop_weight, pop_count_sum = input$pop_weight,
   burned_area = input$fire_weight, burned_pct = input$fire_weight,
   canopypct_mean = 1, balive_sum = 1
  )
  
  X_weighted <- X * variable_weights[colnames(X)]
  return(X_weighted)
 })
 
 # Train SOM reactively
 som_results <- eventReactive(input$run_som, {
  som_grid <- somgrid(xdim = input$xdim, ydim = input$ydim, topo = "hexagonal")
  som_model <- som(as.matrix(weighted_data()), grid = som_grid, rlen = 100)
  
  clusters <- cutree(hclust(dist(som_model$codes[[1]])), k = 5)
  grid$som_cluster <- clusters[som_model$unit.classif]
  
  return(grid)
 })
 
 # Post-hoc analysis: Most common forest type per cluster
 dominant_species <- reactive({
  som_results() %>%
   group_by(som_cluster) %>%
   summarise(
    most_common_spp1 = names(which.max(table(dom_spp1))),
    most_common_spp2 = names(which.max(table(dom_spp2))),
    most_common_spp3 = names(which.max(table(dom_spp3)))
   )
 })
 
 # Render updated SOM map
 output$som_map <- renderPlot({
  ggplot(som_results()) +
   geom_sf(aes(fill = as.factor(som_cluster)), color = NA) +
   scale_fill_viridis_d(name = "SOM Cluster") +
   theme_minimal() +
   labs(title = "Updated SOM Clusters")
 })
 
 # Render dominant species table
 output$dominant_species_table <- renderDT({
  datatable(dominant_species())
 })
 
 # Render dominant species bar plot
 output$dominant_species_plot <- renderPlot({
  dominant_species() %>%
   tidyr::pivot_longer(cols = starts_with("most_common_spp"), names_to = "rank", values_to = "species") %>%
   count(som_cluster, species) %>%
   ggplot(aes(x = as.factor(som_cluster), y = n, fill = species)) +
   geom_bar(stat = "identity", position = "dodge") +
   theme_minimal() +
   labs(title = "Dominant Forest Types by Cluster", x = "SOM Cluster", y = "Count", fill = "Species")
 })
}

library(shiny)

ui <- fluidPage(
 titlePanel("Interactive SOM Clustering for Aspen Fire Management"),
 
 sidebarLayout(
  sidebarPanel(
   sliderInput("wui_weight", "WUI Weight:", min = 1, max = 10, value = 5),
   sliderInput("pop_weight", "Population Weight:", min = 1, max = 10, value = 5),
   sliderInput("fire_weight", "Fire Risk Weight:", min = 1, max = 10, value = 1),
   sliderInput("xdim", "SOM Grid X Dimension:", min = 5, max = 20, value = 10),
   sliderInput("ydim", "SOM Grid Y Dimension:", min = 5, max = 20, value = 10),
   actionButton("run_som", "Update Clustering")
  ),
  
  mainPanel(
   plotOutput("som_map"),
   h3("Dominant Forest Types by Cluster"),
   DTOutput("dominant_species_table"),
   plotOutput("dominant_species_plot")
  )
 )
)