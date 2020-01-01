
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("Simulation Set-up",
                      sidebarLayout(
                        sidebarPanel(htmlOutput("text_vals")),
                        mainPanel(h4("STUFF"),
                                  sliderInput(inputId = "mutationrate", label = "Mutation Rate",
                                              min = 0.005, max = 0.5, value = 0.05),
                                  numericInput(inputId = "maxdivisions", label = "Max Allowed Divisions",
                                               min = 1, max = 10, value = 5, step = 1),
                                  numericInput(inputId = "numrounds", label = "Number of Rounds",
                                               min = 5, max = 100, value = 10, step = 5),
                                  actionButton(inputId = "button1", label = "Submit Values")
                        ))),
             tabPanel("Simulate Single Patient",
                      sidebarLayout(
                        sidebarPanel(h4("STUFF"),
                                     actionButton(inputId = "simulate1", label = "Run Simulation"),
                                     actionButton(inputId = "showmap1", label = "Show Map"),
                                     actionButton(inputId = "button2", label = "Push Button")),
                        mainPanel(plotOutput("PLOT2"))
                      )
             ),
             tabPanel("Simulate Multiple Patient",
                      sidebarLayout(
                        sidebarPanel(h4("STUFF"),
                                     actionButton(inputId = "button3", label = "Push Button")),
                        mainPanel(plotOutput("PLOT3"))
                      )
             )
  )
)
