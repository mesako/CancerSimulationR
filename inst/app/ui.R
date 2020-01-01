
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("PANEL 1",
                      sidebarLayout(
                        sidebarPanel(h4("STUFF"),
                                     sliderInput(inputId = "mutationrate", label = "Mutation Rate",
                                                 min = 0.001, max = 0.5, value = 0.05),
                                     numericInput(inputId = "maxdivisions", label = "Max Allowed Divisions",
                                                 min = 1, max = 10, value = 5, step = 1),
                                     numericInput(inputId = "numrounds", label = "Number of Rounds",
                                                  min = 5, max = 100, value = 10, step = 1),
                                     actionButton(inputId = "button1", label = "Push Button")
                        ),
                        mainPanel(plotOutput("PLOT1")))),
             tabPanel("PANEL 2",
                      sidebarLayout(
                        sidebarPanel(h4("STUFF"),
                                     actionButton(inputId = "button2", label = "Push Button")),
                        mainPanel(plotOutput("PLOT2"))
                      )
             )
  )
)
