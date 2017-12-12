require(shiny)
library(DT)

ui <- shinyUI(
  DT::dataTableOutput('mytable')
)

dat <- data.frame(
  country = c('USA', 'China'),
  flag = c('<img src="american-flag-large.png"></img>',
           '<img src="china.png"></img>'
  )
)

server <- shinyServer(function(input, output){
  output$mytable <- DT::renderDataTable({
    DT::datatable(dat, escape = FALSE)
  })
})

shinyApp(ui, server)
