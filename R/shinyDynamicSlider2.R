#rm(list = ls())
library(shiny)
app <- shinyApp(
  ui = bootstrapPage(
    uiOutput("myList"),
    actionButton("acb1", "Change Value"),
    actionButton("acb2", "Change Value to Range")
  ),
  server = function(input, output, session) {

    slidertype <- reactiveValues()
    slidertype$type <- "default"

    observeEvent(input$acb1,{slidertype$type <- "normal"})
    observeEvent(input$acb2, {slidertype$type <- "range"})

    output$myList <- renderUI({

      if(slidertype$type == "normal"){
        sliderInput("sld1", min = 0, max = 10, label = "y1", value = 2)
      }
      else if(slidertype$type == "range"){
        sliderInput("sld1", min = 0, max = 10, label = "y1", value = c(2,7))
      }
      else{
        sliderInput("sld1", min = 0, max = 10, label = "y1", value = 5)
      }
    })
  })
runApp(app)
