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
      br(),
      conditionalPanel(
        condition = "output.n",
        actionButton("add_btn", "Add Filter")
      ),
      conditionalPanel(
        condition = "output.table",
        actionButton("rm_btn", "Remove Filter"),
        br(),
        br(),
        textOutput("counter2"),
        br()
      ),
      uiOutput("textbox_ui1")
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
server <- function(input, output, session) {

  data <- reactive({getGEO(input$DataID)})

  tab <- reactive({

    x <- data()

    tbl <- switch(input$tableType, "Gene Expression" = extExp(x,long = input$long), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
    tbl
  })

  #START: mess with adding filters, should probably be after submit

  # Track all user inputs
  AllInputs <- reactive({
    reactiveValuesToList(input)
  })

  #END

  observeEvent(input$submit, {
    ##create function to "render" data table before it is displayed to pass the length of the list back to the UI

    # Track the number of input boxes to render
    counter2 <- reactiveValues(n = 0)

    observeEvent(input$add_btn, { #this is called before the inputs are re-rendered, create reactive function to define inputs
      tbl <- tab()
      counter2$n <- counter2$n + 1
      })

    observeEvent(input$rm_btn, {
      tbl <- tab()
      if (counter2$n > 0) counter2$n <- counter2$n - 1
    })

    output$counter2 <- renderPrint(print(counter2$n))

    #used to remove add button when max filters have been reached
    output$n <- reactive({
      tbl <- tab()
      if(counter2$n >= length(colnames(tbl[[2]]))){
        counter2$n <- length(colnames(tbl[[2]]))
        FALSE
      }else{
        TRUE
      }
      #browser()
    })
    outputOptions(output, 'n', suspendWhenHidden = FALSE)

    #MESSING WITH SELECTIZE
    chc <- reactiveValues(ch = sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
    upIn <- reactive({ #find way to correctly define eventExpr
      tbl <- tab()
      n <- counter2$n
      ch <- chc$ch
      browser()
      ch <- match(ch,colnames(tbl[[2]]))
      ch <- ch[!is.na(ch)]
      tbl[[2]][,ch]
      browser()
      lapply(seq_len(n), function(i) {
        updateSelectizeInput(session, inputId = paste0("textin", i), choices = ch)
      })
    })

    textboxes1 <- reactive({ #FIND WHAT MAKES THIS BE CALLED WHEN ADD OR REMOVE BUTTONS ARE HIT
      n <- counter2$n
      tbl <- tab()

      if (n > 0) {
        if(length(tbl[[2]][,1]) <= 50){ #whenever n or tbl are changed, the inputs are re-rendered
          lapply(seq_len(n), function(i) {
            #browser()
            list(
              selectInput(inputId = paste0("selctin", i),
                          label = paste0("Filter by:", i),
                          choices = colnames(tbl[[2]])), #choices are going to be columns of the table, could call tab() now, might have to create diff tab for rendering choices for dataset and here
              selectizeInput(inputId = paste0("textin", i),
                             label = paste0("Value:", i),
                             choices = tbl[[2]][,1]) #cannot have this function call a reactive variable, causes re-render just like changing table type, or add filter
            )
          })
        }else{
          lapply(seq_len(n), function(i) {
            #browser()
            list(
              selectInput(inputId = paste0("selctin", i),
                          label = paste0("Filter by:", i),
                          choices = colnames(tbl[[2]])), #choices are going to be columns of the table, could call tab() now, might have to create diff tab for rendering choices for dataset and here
              textInput(inputId = paste0("textin", i),
                             label = paste0("Value:", i))
            )
          })
        }
      }
    })

    output$textbox_ui1 <- renderUI({ textboxes1() })

    output$length <- renderUI({
        tbl <- tab()
        #browser()
        selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    output$table <- DT::renderDataTable({

      tbl <- tab()
      num <- counter2$n

      #browser()

      if(num > 0){
        y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
        y <- match(y,colnames(tbl[[2]]))
        y <- y[!is.na(y)]
        y <- data.frame(tbl[[2]][,y])
        #look at AllInputs to see if this helps the reset bug

        z <- data.frame(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
        #browser()
        if(identical(z,character(0))){
          z <- 0
        }else{
          z <- row.match(z, y) #could make this more generalized and complicated with grep
        }
      }else{
        z <- 0
      }

      #browser()

      #code to change "page" or dataset
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

      if(num==0 | identical(z,0) | identical(z,integer(0)))  {
        tbl
      }else{
        tbl[z,]
      }
    }, rownames = FALSE)
  })

  observe({
    input$DataID
    counter$countervalue <- 1 #sets to one whenever DataID is changed
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
    y <- character(0)
    z <- character(0)

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

    if(is.na(y)){ #will have to change if statements so that the list of parameters chosen is added to end of file
                  #should simplify if statements, but complicate pasting
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
}

shinyApp(ui = ui, server = server)
