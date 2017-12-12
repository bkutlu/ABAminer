# shiny make slider vale dynamic
library(shiny)
runApp(list(
  ui = bootstrapPage(
    numericInput('n', 'Maximum of slider', 100),
    uiOutput("slider"),
    textOutput("test")
  ),
  server = function(input, output) {
    output$slider <- renderUI({
      sliderInput("myslider", "Slider text", 1,
                  max(input$n, isolate(input$myslider)), 21)
    })

    output$test <- renderText({input$myslider})
  }
))
