require(shiny)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(maps)

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

  observeEvent(input$single_values, {
    output$single_sim_button <- renderUI(
      actionButton(inputId = "single_simulate", label = "Run Simulation")
    )
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

  patient_status <- eventReactive(input$single_simulate, {
    str0 <- "<br/>"
    str1 <- paste0(single_output[[8]], collapse = "<br/>")
    str2 <- single_output[[7]]
    str3 <- "<br/>"
    paste(str0, str1, str2, sep = "<br/>")
  })
  output$patient_status <- renderUI({
    HTML(patient_status())
  })

  show_starting_map <- eventReactive(input$single_simulate, {
    print(starting_map + ggtitle("Starting Cell Map"))
  })
  output$show_starting_map <- renderPlot({
    show_starting_map()
  })

  show_final_map <- eventReactive(input$single_simulate, {
    print(final_map + ggtitle("Final Cell Map"))
  })
  output$show_final_map <- renderPlot({
    show_final_map()
  })

  single_num_plot <- eventReactive(input$single_simulate, {
    print(CancerSimulationR:::PlotTotalCellNum(single_output[[5]]) +
            ggtitle("Number of Cancer Cells in Each Round"))
  })
  output$single_num_plot <- renderPlot({
    single_num_plot()
  })

  single_mut_plot <- eventReactive(input$single_simulate, {
    print(CancerSimulationR:::PlotAverageMutNum(single_output[[6]]) +
            ggtitle("Average Number of Mutations in Each Round"))
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
    multi_text_values$patients <- input$multi_numpatients
    map_dim <<- c(10, 10)
    mutation_set <<- CancerSimulationR:::EstablishMutations()
    patient_preaction_set <<- CancerSimulationR:::EstablishPatientPreActions()
    patient_postaction_set <<- CancerSimulationR:::EstablishPatientPostActions()
    mutation_encoding <<- CancerSimulationR:::MakeMutationEncoding(mutation_set)
  })
  output$show_multi_text <- renderUI({
    str0 <- "<br/>"
    str1 <- paste("Mutation rate assigned to", multi_text_values$mutrate)
    str2 <- paste("Max number of cell divisions to", multi_text_values$maxdiv)
    str3 <- paste("Number of rounds assigned to", multi_text_values$rounds)
    str4 <- paste("Number of patients to simulate assigned to", multi_text_values$patients)
    str5 <- "<br/>"
    HTML(paste(str0, str1, str2, str3, str4, str5, sep = "<br/>"))
  })

  observeEvent(input$multi_values, {
    output$multi_sim_button <- renderUI(
      actionButton(inputId = "multi_simulate", label = "Run Simulations")
    )
  })

  ### RUN MULTIPLE SIMULATION
  observeEvent(input$multi_simulate, {
    multi_output <<- list()
    for (i in 1:multi_text_values$patients) {
      starting_map <<- CancerSimulationR:::GenerateMap(map_dim)
      cell_mut_rate <<- CancerSimulationR:::GenerateMutationMatrix(map_dim, starting_map, multi_text_values$mutrate)
      cell_divisions <<- CancerSimulationR:::GenerateDivisionMatrix(map_dim, starting_map)
      cell_states <<- CancerSimulationR:::MakeStartingCellStateMatrix(map_dim)
      multi_output[[i]] <<- CancerSimulationR:::RunPatientSimulation(multi_text_values$rounds, starting_map, mutation_set,
                                                                     patient_preaction_set, patient_postaction_set,
                                                                     mutation_encoding, cell_mut_rate, map_dim,
                                                                     cell_divisions, cell_states, multi_text_values$maxdiv)

    }
  })

  show_survival_curve <- eventReactive(input$multi_simulate, {
    patient.summary.statements <<- c()
    for (i in 1:multi_text_values$patients) {
      patient.summary.statements <<- c(patient.summary.statements,
                                       multi_output[[i]][[7]])
    }
    print(CancerSimulationR:::PlotMultiSurvivalCurve(patient.summary.statements,
                                                     multi_text_values$rounds) +
            ggtitle("Patient Survival Curve"))
  })
  output$show_survival_curve <- renderPlot({
    show_survival_curve()
  })


  show_patient_status <- eventReactive(input$multi_simulate, {
    patient.summary.statements <<- c()
    for (i in 1:multi_text_values$patients) {
      patient.summary.statements <<- c(patient.summary.statements,
                                       multi_output[[i]][[7]])
    }
    print(CancerSimulationR:::PlotMultiPatientSummary(patient.summary.statements) +
            ggtitle("Number of Patients Grouped by Status"))
  })
  output$show_patient_status <- renderPlot({
    show_patient_status()
  })

  multi_mut_plot <- eventReactive(input$multi_simulate, {
    all.average.mut.nums <<- c()
    for (i in 1:multi_text_values$patients) {
      temp <- multi_output[[i]][[6]]
      if (length(temp) < multi_text_values$rounds) {
        temp <- c(temp, rep(temp[length(temp)],
                            times = (multi_text_values$rounds - length(temp))))
      }
      all.average.mut.nums <<- rbind(all.average.mut.nums, temp)
    }
    print(CancerSimulationR:::MultiPlotAverageMutNum(all.average.mut.nums)+
            ggtitle("Average Number of Mutations in Each Round"))
  })
  output$multi_mut_plot <- renderPlot({
    multi_mut_plot()
  })

  multi_num_plot <- eventReactive(input$multi_simulate, {
    all.total.cell.nums <<- c()
    for (i in 1:multi_text_values$patients) {
      temp <- multi_output[[i]][[5]]
      if (length(temp) < multi_text_values$rounds) {
        temp <- c(temp, rep(temp[length(temp)],
                            times = (multi_text_values$rounds - length(temp))))
      }
      all.total.cell.nums <<- rbind(all.total.cell.nums, temp)
    }
    print(CancerSimulationR:::MultiPlotTotalCellNum(all.total.cell.nums) +
            ggtitle("Number of Cancer Cells in Each Round"))
  })
  output$multi_num_plot <- renderPlot({
    multi_num_plot()
  })
}
