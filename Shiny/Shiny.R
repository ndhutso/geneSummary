library(shiny)
library(geneSummary)

# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
ui <- fluidPage(

  # App title ----
  titlePanel("Data Tables"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      h1("Select data and type of table"),
      selectInput(inputId = "tableType", label = "Type of data table", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations")),
      textInput(inputId = "DataID", label = "GEO accession ID:")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Data table ----
      tableOutput("table")

    )
  )
)
# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot

  table <- reactive({

    x <- getGEO(input$DataID)
    switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x) )

  })

  output$table <- renderTable({

    table()

  })

}

shinyApp(ui = ui, server = server)
