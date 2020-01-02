
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("Simulate Single Patient",
                      sidebarLayout(
                        sidebarPanel(h4("Set-up Single Simulation"),
                                     sliderInput(inputId = "single_mutationrate", label = "Mutation Rate",
                                                 min = 0.005, max = 0.5, value = 0.05),
                                     numericInput(inputId = "single_maxdivisions", label = "Max Allowed Divisions",
                                                  min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "single_numrounds", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 5),
                                     actionButton(inputId = "single_values", label = "Submit Values"),
                                     htmlOutput("show_single_text"),
                                     actionButton(inputId = "single_simulate", label = "Run Simulation")),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Patient Summary",
                                     htmlOutput("patient_status")),
                            tabPanel("Cell Maps", splitLayout(cellWidths = c("50%", "50%"),
                                                              plotOutput("show_starting_map"),
                                                              plotOutput("show_final_map"))),
                            tabPanel("Summary Plots", splitLayout(cellWidths = c("50%", "50%"),
                                                                  plotOutput("single_num_plot"),
                                                                  plotOutput("single_mut_plot"))),
                            tabPanel("Concept Questions",
                                     h3("QUESTIONS"),
                                     p("1. first question"),
                                     p("2. second question"),
                                     p("3. third question")
                            )
                          )
                        )
                      )
             ),
             tabPanel("Simulate Multiple Patient",
                      sidebarLayout(
                        sidebarPanel(h4("Set-up Multiple Simulations"),
                                     sliderInput(inputId = "multi_mutationrate", label = "Mutation Rate",
                                                 min = 0.005, max = 0.5, value = 0.05),
                                     numericInput(inputId = "multi_maxdivisions", label = "Max Allowed Divisions",
                                                  min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "multi_numrounds", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 5),
                                     actionButton(inputId = "multi_values", label = "Submit Values"),
                                     htmlOutput("show_multi_text"),
                                     actionButton(inputId = "multi_simulate", label = "Run Simulations")),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Patient Summaries",
                                     h3("INFO")),
                            tabPanel("Example Cell Maps", splitLayout(cellWidths = c("50%", "50%"),
                                                                      h3("BOOP"),
                                                                      h3("BOOP"))),
                            tabPanel("Summary Plots", splitLayout(cellWidths = c("50%", "50%"),
                                                                  h3("BOOP"),
                                                                  h3("BOOP"))),
                            tabPanel("Concept Questions",
                                     h3("QUESTIONS"),
                                     p("1. first question"),
                                     p("2. second question"),
                                     p("3. third question")
                            )
                          )
                        )

                      )
             )
  )
)
