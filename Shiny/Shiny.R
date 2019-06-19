library(DT)
library(shiny)
library(geneSummary)

# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
ui <- fluidPage(

  titlePanel("Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations"), selected = "Gene Expression"),
      textInput(inputId = "DataID", label = "GEO accession ID:"),
      textInput(inputId = "DataName", label = "Data set name:", placeholder = "All"),
      conditionalPanel(
        condition = "input.tableType != 'Sample Annotations'",
        textInput(inputId = "geneSymbol", label = "Gene Symbol:", placeholder = "All")
      ),
      actionButton("submit", label = "Submit"),
      br(),
      br(),
      strong("Graphs for GSE43452:"),
      br(),
      actionButton("graph1", label = "TP53 Concentration"),
      br(),
      br(),
      actionButton("graph2", label = "TP53 in DBTRG and U87"),
      br(),
      br(),
      actionButton("graph3", label = "Top Variance Expressions"),
      br(),
      br(),
      actionButton("reset", label = "Reset")
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      conditionalPanel(
        condition = "output.table",
        conditionalPanel(
          condition = "input.tableType == 'Gene Expression'",
          titlePanel("Gene Expression Data")
        ),
        conditionalPanel(
          condition = "input.tableType == 'Gene Annotations'",
          titlePanel("Gene Annotations")
        ),
        conditionalPanel(
          condition = "input.tableType == 'Sample Annotations'",
          titlePanel("Sample Annotations")
        ),
        #page selector, have it default to 1 so the first element of the list will be ouput by render graph at first
        selectInput(inputId = "page", label = "Page:",choices = 1:10, selected = 1)
      ),

      #Hide errors
      #tags$style(type="text/css",
      #           ".shiny-output-error { visibility: hidden; }",
      #           ".shiny-output-error:before { visibility: hidden; }"),
      # Output: Data table ----
      DT::dataTableOutput("table"),

      conditionalPanel(
        condition = "output.table",
        actionButton("save", label = "Save Data")
      ),

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

  observeEvent(input$submit, {
    output$table <- DT::renderDataTable({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      z <- strsplit(input$DataName,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
        if(identical(z, character(0))){
          t <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
        }else{
          t <- switch(input$tableType, "Gene Expression" = extExp(x,dName = z), "Gene Annotations" = extGene(x,dName = z),"Sample Annotations" = extSample(x,z))
        }
      }else if(identical(z, character(0))){
        t <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }else{
        t <- switch(input$tableType, "Gene Expression" = extExp(x,y,z), "Gene Annotations" = extGene(x,y,z),"Sample Annotations" = extSample(x,z))
      }

      if(length(t)>1){
        n <- as.numeric(input$page)
        t <- t[[n]]
      }
      t

      ##PROBLEM: GEOQuery doesnt like to be imported anymore
    }, rownames = TRUE)
  })

  observeEvent(input$save, {
    x <- data()

    y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

    z <- strsplit(input$DataName,", ",fixed = TRUE)[[1]]

    D1a <- extGene(x,y,z)
    D2a <- extExp(x,y,z)
    D2b <- extSample(x,z)

    dir.create("data")
    setwd("data")
    save(D1a, file = "geneAnnotation.RData")
    save(D2a, file = "geneExpression.RData")
    save(D2b, file = "sampleAnnotation.RData")

  })

  ##ONLY WORKS FOR GSE43452
  #use reactive values to define graph and just call the variable in render plot
  v <- reactiveValues(data = NULL)
  dataDef <- reactive({getGEO("GSE43452")})

  observeEvent(input$graph1, {
    x <- dataDef()

    y <- NA

    z <- NA

    D1a <- extGene(x,y,z)
    D2a <- extExp(x,y,z)
    D2b <- extSample(x,z)

    v$graph <- bar(D1a,D2a,D2b)
  })
  observeEvent(input$graph2, {
    x <- dataDef()

    y <- NA

    z <- NA

    D1a <- extGene(x,y,z)
    D2a <- extExp(x,y,z)
    D2b <- extSample(x,z)

    v$graph <- box(D1a,D2a)
  })
  observeEvent(input$graph3, {
    x <- dataDef()

    y <- NA

    z <- NA

    D1a <- extGene(x,y,z)
    D2a <- extExp(x,y,z)
    D2b <- extSample(x,z)

    v$graph <- hist(D1a,D2a)
  })

  observeEvent(input$reset, {
    v$graph <- NULL
  })

  output$graph <- renderPlot({

    v$graph

  })

}

shinyApp(ui = ui, server = server)
