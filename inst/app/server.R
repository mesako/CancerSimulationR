library(shiny)
library(RColorBrewer)
library(ggplot2)
library(reshape)

server <- function(input, output) {
  ### SIMULATE SINGLE PATIENT
  ### SIMULATION SET-UP
  single_text_values <- reactiveValues()
  observeEvent(input$single_values, {
    single_text_values$mutrate <- input$single_mutationrate
    single_text_values$maxdiv <- input$single_maxdivisions
    single_text_values$rounds <- input$single_numrounds
    map_dim <<- c(10, 10)
    mutation_set <<- CancerSimulationR:::EstablishMutations()
    patient_preaction_set <<- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <<- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <<- CancerSimulationR:::MakeMutationEncoding(mutation_set)
  })
  output$show_single_text <- renderUI({
    str0 <- "<br/>"
    str1 <- paste("Mutation rate assigned to", single_text_values$mutrate)
    str2 <- paste("Max number of cell divisions to", single_text_values$maxdiv)
    str3 <- paste("Number of rounds assigned to", single_text_values$rounds)
    str4 <- "<br/>"
    HTML(paste(str0, str1, str2, str3, str4, sep = "<br/>"))
  })

  ### RUN SINGLE SIMULATION
  observeEvent(input$single_simulate, {
    starting_map <<- CancerSimulationR:::GenerateMap(map_dim)
    cell_mut_rate <<- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, single_text_values$mutrate)
    cell_divisions <<- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
    cell_states <<- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
    single_output <<- CancerSimulationR:::RunPatientSimulation(single_text_values$rounds, starting_map, mutation_set,
                                                             patient_preaction_set, patient_postaction_set,
                                                             mutation_encoding, cell_mut_rate, map_dim,
                                                             cell_divisions, cell_states, single_text_values$maxdiv)
    final_map <<- single_output[[1]]
  })
  show_starting_map <- eventReactive(input$single_simulate, {
    print(starting_map)
  })
  output$show_starting_map <- renderPlot({
    show_starting_map()
  })

  show_final_map <- eventReactive(input$single_simulate, {
    print(final_map)
  })
  output$show_final_map <- renderPlot({
    show_final_map()
  })

  single_num_plot <- eventReactive(input$single_simulate, {
    print(CancerSimulationR:::PlotTotalCellNum(single_output[[5]]))
  })
  output$single_num_plot <- renderPlot({
    single_num_plot()
  })

  single_mut_plot <- eventReactive(input$single_simulate, {
    print(CancerSimulationR:::PlotAverageMutNum(single_output[[6]]))
  })
  output$single_mut_plot <- renderPlot({
    single_mut_plot()
  })

  ### SIMULATE MULTIPLE PATIENTS
  ### SIMULATION SET-UP
  multi_text_values <- reactiveValues()
  observeEvent(input$multi_values, {
    multi_text_values$mutrate <- input$multi_mutationrate
    multi_text_values$maxdiv <- input$multi_maxdivisions
    multi_text_values$rounds <- input$multi_numrounds
    map_dim <<- c(10, 10)
    starting_map <<- CancerSimulationR:::GenerateMap(map_dim)
    mutation_set <<- CancerSimulationR:::EstablishMutations()
    patient_preaction_set <<- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <<- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <<- CancerSimulationR:::MakeMutationEncoding(mutation_set)
    cell_mut_rate <<- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, multi_text_values$mutrate)
    cell_divisions <<- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
    cell_states <<- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
  })
  output$show_multi_text <- renderUI({
    str0 <- "<br/>"
    str1 <- paste("Mutation rate assigned to", multi_text_values$mutrate)
    str2 <- paste("Max number of cell divisions to", multi_text_values$maxdiv)
    str3 <- paste("Number of rounds assigned to", multi_text_values$rounds)
    str4 <- "<br/>"
    HTML(paste(str0, str1, str2, str3, str4, sep = "<br/>"))
  })

  ### RUN MULTIPLE SIMULATION
  observeEvent(input$multi_simulate, {
    multi_output <<- CancerSimulationR:::RunPatientSimulation(single_text_values$rounds, starting_map, mutation_set,
                                                             patient_preaction_set, patient_postaction_set,
                                                             mutation_encoding, cell_mut_rate, map_dim,
                                                             cell_divisions, cell_states, single_text_values$maxdiv)
  })
  final_map <- eventReactive(input$multi_simulate, {
    print(multi_output[[1]])
  })
  output$final_map <- renderPlot({
    final_map()
  })
}
