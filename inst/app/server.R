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
    map_dim <- c(10, 10)
    print(map_dim)
    starting_map <- CancerSimulationR:::GenerateMap(map_dim)
    print(starting_map)
    mutation_set <- CancerSimulationR:::EstablishMutations()
    print(text_values$mutrate)
    print(text_values$maxdiv)
    print(text_values$rounds)
    patient_preaction_set <- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <- CancerSimulationR:::MakeMutationEncoding(mutation_set)
    cell_mut_rate <- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, text_values$mutrate)
    cell_divisions <- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
    cell_states <- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
  })
  output$text_vals <- renderUI({
    str1 <- paste("Mutation rate assigned to", text_values$mutrate)
    str2 <- paste("Max number of cell divisions to", text_values$maxdiv)
    str3 <- paste("Number of rounds assigned to", text_values$rounds)
    HTML(paste(str1, str2, str3, sep = "<br/>"))
  })

  ### SIMULATE SINGLE PATIENT
  observeEvent(input$simulate1, {
    save_output <- CancerSimulationR:::RunPatientSimulation(text_values$rounds, starting_map, mutation_set,
                                                                    patient_preaction_set, patient_postaction_set,
                                                                    mutation_encoding, cell_mut_rate,
                                                                    cell_divisions, cell_states, text_values$maxdiv)
  })

  PLOT1 <- eventReactive(input$showmap1, {
    print(save_output[[1]])
  })
  output$PLOT1 <- renderPlot({
    PLOT1()
  })


  PLOT2 <- eventReactive(input$button2, {
    # fill in plot
  })
  output$PLOT2 <- renderPlot({
    PLOT2()
  })
}
