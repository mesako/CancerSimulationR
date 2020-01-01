library(shiny)
library(RColorBrewer)
library(ggplot2)
library(reshape)

server <- function(input, output) {
  
  react.vals <- reactiveValues()

  PLOT1 <- eventReactive(input$button1, {
    # fill in plot
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
