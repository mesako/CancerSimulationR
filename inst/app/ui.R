
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("Simulate Single Patient",
                      sidebarLayout(
                        sidebarPanel(h4("Set-up Single Simulation"),
                                     sliderInput(inputId = "mutationrate", label = "Mutation Rate",
                                                 min = 0.005, max = 0.5, value = 0.05),
                                     numericInput(inputId = "maxdivisions", label = "Max Allowed Divisions",
                                                  min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "numrounds", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 5),
                                     actionButton(inputId = "button1", label = "Submit Values"),
                                     htmlOutput("text_vals"),
                                     actionButton(inputId = "simulate1", label = "Run Simulation")),
                        mainPanel(plotOutput("final_map"))
                      )
             ),
             tabPanel("Simulate Multiple Patient",
                      sidebarLayout(
                        sidebarPanel(h4("Set-up Multiple Simulations"),
                                     sliderInput(inputId = "mutationrate2", label = "Mutation Rate",
                                                 min = 0.005, max = 0.5, value = 0.05),
                                     numericInput(inputId = "maxdivisions2", label = "Max Allowed Divisions",
                                                  min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "numrounds2", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 5),
                                     actionButton(inputId = "button2", label = "Submit Values"),
                                     htmlOutput("text_vals2"),
                                     actionButton(inputId = "simulate2", label = "Run Simulations")),
                        mainPanel(plotOutput("summary_plots"))
                      )
             )
  )
)
