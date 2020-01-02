library(shiny)
library(RColorBrewer)
library(ggplot2)
library(reshape)

server <- function(input, output) {
  ### SIMULATION SET-UP
  text_values <- reactiveValues()
  observeEvent(input$button1, {
    text_values$mutrate <- input$mutationrate
    text_values$maxdiv <- input$maxdivisions
    text_values$rounds <- input$numrounds
    map_dim <<- c(10, 10)
    starting_map <<- CancerSimulationR:::GenerateMap(map_dim)
    mutation_set <<- CancerSimulationR:::EstablishMutations()
    patient_preaction_set <<- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <<- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <<- CancerSimulationR:::MakeMutationEncoding(mutation_set)
    cell_mut_rate <<- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, text_values$mutrate)
    cell_divisions <<- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
    cell_states <<- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
  })
  output$text_vals <- renderUI({
    str0 <- "<br/>"
    str1 <- paste("Mutation rate assigned to", text_values$mutrate)
    str2 <- paste("Max number of cell divisions to", text_values$maxdiv)
    str3 <- paste("Number of rounds assigned to", text_values$rounds)
    str4 <- "<br/>"
    HTML(paste(str0, str1, str2, str3, str4, sep = "<br/>"))
  })

  ### SIMULATE SINGLE PATIENT
  observeEvent(input$simulate1, {
    save_output <<- CancerSimulationR:::RunPatientSimulation(text_values$rounds, starting_map, mutation_set,
                                                             patient_preaction_set, patient_postaction_set,
                                                             mutation_encoding, cell_mut_rate, map_dim,
                                                             cell_divisions, cell_states, text_values$maxdiv)
  })
  final_map <- eventReactive(input$simulate1, {
    print(save_output[[1]])
  })
  output$final_map <- renderPlot({
    final_map()
  })


  ### SIMULATE MULTIPLE PATIENTS
  text_values2 <- reactiveValues()
  observeEvent(input$button2, {
    text_values2$mutrate <- input$mutationrate2
    text_values2$maxdiv <- input$maxdivisions2
    text_values2$rounds <- input$numrounds2
    map_dim <<- c(10, 10)
    starting_map <<- CancerSimulationR:::GenerateMap(map_dim)
    mutation_set <<- CancerSimulationR:::EstablishMutations()
    patient_preaction_set <<- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <<- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <<- CancerSimulationR:::MakeMutationEncoding(mutation_set)
    cell_mut_rate <<- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, text_values$mutrate)
    cell_divisions <<- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
    cell_states <<- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
  })
  output$text_vals2 <- renderUI({
    str0 <- "<br/>"
    str1 <- paste("Mutation rate assigned to", text_values2$mutrate)
    str2 <- paste("Max number of cell divisions to", text_values2$maxdiv)
    str3 <- paste("Number of rounds assigned to", text_values2$rounds)
    str4 <- "<br/>"
    HTML(paste(str0, str1, str2, str3, str4, sep = "<br/>"))
  })
}
