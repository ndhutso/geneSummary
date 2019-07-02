# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
#GEOquery, stringr
ui <- fluidPage(

  titlePanel("Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations"), selected = "Gene Expression"),
      conditionalPanel(
        condition = "input.tableType == 'Gene Expression'",
        checkboxInput("long", "Long format data", value = FALSE)
      ),
      textInput(inputId = "DataID", label = "GEO accession ID:"),
      actionButton("submit", label = "Submit"),
      br(),
      uiOutput("textbox_ui"),
      br(),
      actionButton("add_btn", "Add Filter"),
      actionButton("rm_btn", "Remove Filter"),
      br(),
      textOutput("counter"),
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
        )
      ),

      uiOutput("length"),

      conditionalPanel(
        condition = "output.table",
        actionButton("save", label = "Save Data")
      ),
      br(),

      #Hide errors
      #tags$style(type="text/css",
      #           ".shiny-output-error { visibility: hidden; }",
      #           ".shiny-output-error:before { visibility: hidden; }"),
      # Output: Data table ----
      DT::dataTableOutput("table"),

      plotOutput("graph")

    )
  )
)
# Define server logic
server <- function(input, output) {

  #START: ess with adding filters



  #END


  data <- reactive({getGEO(input$DataID)})

  observeEvent(input$submit, {
    ##create function to "render" data table before it is displayed to pass the length of the list back to the UI
    output$length <- renderUI({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      z <- strsplit(input$sampleName,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
        if(identical(z, character(0))){
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
        }else{
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x,z))
        }
      }else if(identical(z, character(0))){
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x,z))
      }

        #browser()
        selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    output$table <- DT::renderDataTable({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      z <- strsplit(input$sampleName,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
        if(identical(z, character(0))){
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x,long = input$long), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
        }else{
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x,long = input$long), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x,z))
        }
      }else if(identical(z, character(0))){
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y,long = input$long), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y,long = input$long), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x,z))
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

    }, rownames = FALSE)
  })

  observe({
    input$DataID
    counter$countervalue <- 1 #sets to one whenever DataID is changed
    #setwd(counter$choice) #have to pass choice to here
  })

  counter <- reactiveValues(countervalue = 0,choice = getwd())

  observeEvent(input$save, {

    #browser()

    if(counter$countervalue == 1){
      if (exists('utils::choose.dir')) {
        counter$choice <- choose.dir(caption = "Select folder")
      } else {
        counter$choice <- tk_choose.dir(caption = "Select folder")
      }
      setwd(counter$choice)
        #only call choose.dir if input$DataID is changed after save is hit

      d <- paste("data.",input$page,sep = "")
      dir.create(d)
      setwd(d)
    }else if(counter$countervalue > 1){
      setwd(counter$choice)
      d <- paste("data.",input$page,sep = "")
      dir.create(d)
      setwd(d)
    }

    counter$countervalue <- counter$countervalue + 1

    #browser()

    x <- data()

    y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

    z <- strsplit(input$sampleName,", ",fixed = TRUE)[[1]]

    if(identical(y, character(0))){
      if(identical(z, character(0))){
        y <- NA
        z <- NA
      }else{
        y <- NA
      }
    }else if(identical(z, character(0))){
      z <- NA
    }

    if(is.na(y)){
      if(is.na(z)){
        if(input$tableType == "Gene Expression"){
          tbl <- extExp(x,y,long = input$long)
          #browser()
          if(input$long){
            saveRDS(tbl, file = paste("geneExpression.long.rds", sep=""))
          }else{
            saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
          }
        }else if(input$tableType == "Gene Annotations"){
          tbl <- extGene(x,y)
          saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
        }else{
          #browser()
          tbl <- extSample(x,z)
          saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
        }
      }else{
        if(input$tableType == "Gene Expression"){
          tbl <- extExp(x,y,long = input$long)
          #browser()
          if(input$long){
            saveRDS(tbl, file = paste("geneExpression.long.rds", sep=""))
          }else{
            saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
          }
        }else if(input$tableType == "Gene Annotations"){
          tbl <- extGene(x,y)
          saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
        }else{
          #browser()
          tbl <- extSample(x,z)
          saveRDS(tbl, file = paste("sampleAnnotation.",paste(z,collapse = "."),".rds", sep=""))
        }
      }
    }else if(is.na(z)){
      if(input$tableType == "Gene Expression"){
        tbl <- extExp(x,y,long = input$long)
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".long.rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        tbl <- extGene(x,y)
        saveRDS(tbl, file = paste("geneAnnotation.",paste(y,collapse = "."),".rds", sep=""))
      }else{
        #browser()
        tbl <- extSample(x,z)
        saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
      }
    }else{
      if(input$tableType == "Gene Expression"){
        tbl <- extExp(x,y,long = input$long)
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".long.rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        tbl <- extGene(x,y)
        saveRDS(tbl, file = paste("geneAnnotation.",paste(y,collapse = "."),".rds", sep=""))
      }else{
        #browser()
        tbl <- extSample(x,z)
        saveRDS(tbl, file = paste("sampleAnnotation.",paste(z,collapse = "."),".rds", sep=""))
      }
    }
  })

  ##ONLY WORKS FOR GSE43452
  #use reactive values to define graph and just call the variable in render plot
  v <- reactiveValues(data = NULL)
  dataDef <- reactive({getGEO("GSE43452")})

  observeEvent(input$graph1, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
    D2b <- extSample(x)

    v$graph <- bar(D1a[[2]],D2a[[2]],D2b[[2]])
  })
  observeEvent(input$graph2, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
    D2b <- extSample(x)

    v$graph <- box(D1a[[2]],D2a[[2]],D2b[[2]])
  })
  observeEvent(input$graph3, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
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
