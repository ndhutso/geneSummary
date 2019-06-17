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

      h3("Select data and type of table"),
      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations")),
      textInput(inputId = "geneSymbol", label = "Gene Symbol:", placeholder = "All"),
      textInput(inputId = "DataName", label = "Data set name:", placeholder = "All"),
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
  data <- reactive({getGEO(input$DataID)})
  symbol <- reactive({strsplit(input$geneSymbol,", ")[[1]]})
  name <- reactive({strsplit(input$DataName,", ")[[1]]})

  table <- reactive({

    x <- data()
    y <- symbol()
    z <- name()

    switch(input$tableType, "Gene Expression" = extExp(x,y,z), "Gene Annotations" = extGene(x,y,z),"Sample Annotations" = extSample(x,z) )

  })

  output$table <- renderTable({

    table()

  })

}

shinyApp(ui = ui, server = server)
