library(shiny)
library(geneSummary)

# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
ui <- fluidPage(

  # App title ----
  titlePanel("Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations"), selected = "Gene Expression"),
      textInput(inputId = "DataID", label = "GEO accession ID:"),
      textInput(inputId = "DataName", label = "Data set name:", placeholder = "All"),
      textInput(inputId = "geneSymbol", label = "Gene Symbol:", placeholder = "All"),
      strong("Graph:"),
      br(),
      actionButton("graph1", label = "TP53 Concentration"),
      br(),
      br(),
      actionButton("graph2", label = "TP53 in DBTRG and U87"),
      br(),
      br(),
      actionButton("graph3", label = "Top Variance Expressions")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      # Output: Data table ----
      tableOutput("table"),
      plotOutput("graph")

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
  symbol <- reactive({strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]})
  name <- reactive({strsplit(input$DataName,", ",fixed = TRUE)[[1]]})

  table <- reactive({

    x <- data()
    y <- symbol()
    z <- name()

    switch(input$tableType, "Gene Expression" = extExp(x,y,z), "Gene Annotations" = extGene(x,y,z),"Sample Annotations" = extSample(x,z))
    ##PROBLEM:  DOESNT LIKE THE IF AND ELSE STATEMENTS
  })

  output$table <- renderTable({

    table()

  })

  output$graph <- renderPlot({

    switch(input$action)

  })

}

shinyApp(ui = ui, server = server)
