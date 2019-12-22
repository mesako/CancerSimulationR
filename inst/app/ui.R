
ui <- fluidPage(
  navbarPage("Cancer Simulation",
             tabPanel("PANEL 1",
                      sidebarLayout(
                        sidebarPanel(h4("STUFF"),
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
