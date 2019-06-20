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

      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),

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
        )#,
        #page selector, have it default to 1 so the first element of the list will be ouput by render graph at first

      ),

      uiOutput("length"),

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
    ##create function to "render" data table before it is displayed to pass the length of the list back to the UI
    output$length <- renderUI({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }

        #browser()
        selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    output$table <- DT::renderDataTable({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }

      #browser()
      if(!is.data.frame(tbl[[2]])){
          n <- input$page #being delayed so is passing a NULL
          n <- match(n, tbl[[1]])
          #browser()
          if(is.null(n)){
            n <- 1
          }
          #browser()
          tbl <- tbl[[2]][[n]]
      }else{
        tbl <- tbl[[2]]
      }

      tbl

    }, rownames = TRUE)
  })

  observeEvent(input$save, {
    x <- data()
    y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

    if(identical(y, character(0))){
        y <- NA
    }

    D1a <- extGene(x,y)
    D2a <- extExp(x,y)
    D2b <- extSample(x)

    browser()

    dir.create("data")
    setwd("data")
    saveRDS(D1a, file = "geneAnnotation.rds")
    saveRDS(D2a, file = "geneExpression.rds")
    saveRDS(D2b, file = "sampleAnnotation.rds")

  })

  ##ONLY WORKS FOR GSE43452
  #use reactive values to define graph and just call the variable in render plot
  v <- reactiveValues(data = NULL)
  dataDef <- reactive({getGEO("GSE43452")})

  observeEvent(input$graph1, {
    x <- dataDef()

    y <- NA

    D1a <- extGene(x,y)
    D2a <- extExp(x,y)
    D2b <- extSample(x)

    v$graph <- bar(D1a[[2]],D2a[[2]],D2b[[2]])
  })
  observeEvent(input$graph2, {
    x <- dataDef()

    y <- NA

    D1a <- extGene(x,y)
    D2a <- extExp(x,y)
    D2b <- extSample(x)

    v$graph <- box(D1a[[2]],D2a[[2]])
  })
  observeEvent(input$graph3, {
    x <- dataDef()

    y <- NA

    D1a <- extGene(x,y)
    D2a <- extExp(x,y)
    D2b <- extSample(x)

    v$graph <- hist(D1a[[2]],D2a[[2]])
  })

  observeEvent(input$reset, {
    v$graph <- NULL
  })

  output$graph <- renderPlot({

    v$graph

  })

}

shinyApp(ui = ui, server = server)
