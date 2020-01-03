
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("Simulate Single Patient",
                      sidebarLayout(
                        sidebarPanel(h4("Set-up Single Simulation"),
                                     sliderInput(inputId = "single_mutationrate", label = "Mutation Rate",
                                                 min = 0.005, max = 0.5, value = 0.05),
                                     numericInput(inputId = "single_maxdivisions", label = "Max Allowed Cell Divisions",
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
                                     # h3("QUESTIONS"),
                                     br(), br(),
                                     p("1. How does changing the different parameters affect the final patient status?"),
                                     p("2. How does changing the different parameters affect the time it takes for the patient to die (if at all)?"),
                                     p("3. How does changing the different parameters affect the number of cancer cells across all of the rounds?"),
                                     p("4. How dodoes changing the different parameters affect the average number of mutations across all of the rounds?")
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
                                     numericInput(inputId = "multi_maxdivisions", label = "Max Allowed Cell Divisions",
                                                  min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "multi_numrounds", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 5),
                                     numericInput(inputId = "multi_numpatients", label = "Number of Patients",
                                                  min = 2, max = 20, value = 10, step = 2),
                                     actionButton(inputId = "multi_values", label = "Submit Values"),
                                     htmlOutput("show_multi_text"),
                                     actionButton(inputId = "multi_simulate", label = "Run Simulations")),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Patient Summaries", verticalLayout(
                              plotOutput("show_patient_status"),
                              plotOutput("show_survival_curve"))),
                            tabPanel("Summary Plots", verticalLayout(
                              plotOutput("multi_num_plot"),
                              plotOutput("multi_mut_plot"))),
                            tabPanel("Concept Questions",
                                     # h3("QUESTIONS"),
                                     br(), br(),
                                     p("1. How do the different parameters affect the distribution of final patient statuses?"),
                                     p("2. How do the different parameters affect the time it takes for patients to die (if at all)?"),
                                     p("3. How does the number of cancer cells compare between patients, and across all of the rounds?"),
                                     p("4. How does the average number of mutations compare between patients, and across all of the rounds?")
                            )
                          )
                        )

                      )
             )
  )
)
