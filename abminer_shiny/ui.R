library(shiny)
library(DT)

# Define UI for dataset viewer application
fluidPage(

  # Application title
  titlePanel("Allen Brain Atlas Hypothalamus Specific Genes"),

  # Sidebar with controls to provide a caption, select a dataset,
  # and specify the number of observations to view. Note that
  # changes made to the caption in the textInput control are
  # updated in the output area immediately as you type
  sidebarLayout(
    sidebarPanel(

      selectInput("structure", "Choose a Hypothalamus region:",
                  choices = c("VMH", "PVHD", "ARC", "LHA"))

      # numericInput("obs", "Number of observations to view:", 10)
    ),


    # Show the caption, a summary of the dataset and an HTML
    # table with the requested number of observations
    mainPanel(

      # tableOutput("view")
     DT::dataTableOutput('view')
    )
  )
)
