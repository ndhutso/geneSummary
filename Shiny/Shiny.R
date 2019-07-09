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
        br()
      ),
      uiOutput("textbox_ui1") #displays the filter list output
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      #prevents background warnings from being output for the user to see
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

      # Output: Data table ----
      DT::dataTableOutput("table"),

    )
  )
)
# Define server logic
server <- function(input, output, session) {

  data <- reactive({getGEO(input$DataID)})

  tab <- reactive({ #pre-rendered data and table so that the different types of tables can be accessed quickly without having to re-load the data from online

    x <- data()

    tbl <- switch(input$tableType, "Gene Expression" = extExp(x,long = input$long), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
    tbl
  })

  # Track all user inputs so that they are not reset when the number of filters is changed (they used to get re-rendered blank)
  AllInputs <- reactive({
    reactiveValuesToList(input)
  })

  observeEvent(input$submit, {

    #track the number of input boxes to render
    counter2 <- reactiveValues(n = 0)

    observeEvent(input$add_btn, {
      counter2$n <- counter2$n + 1
      })

    observeEvent(input$rm_btn, {
      if (counter2$n > 0) counter2$n <- counter2$n - 1
    })

    output$counter2 <- renderPrint(print(counter2$n))

    #used to remove add button when max filters have been reached
    output$n <- reactive({
      tbl <- tab()
      if(input$tableType != "Gene Expression"){
        if(counter2$n >= length(colnames(tbl[[2]]))){
          counter2$n <- length(colnames(tbl[[2]]))
          FALSE
        }else{
          TRUE
        }
      }else{ #REMOVE WHEN NUMERIC FILTERS IMPLEMENTED
        if(counter2$n>=2){
          counter2$n <- 2
          FALSE
        }else{
          TRUE
        }
      }

      #browser()
    })
    outputOptions(output, 'n', suspendWhenHidden = FALSE)

    #create UI inputs
    textboxes1 <- reactive({
      n <- counter2$n
      tbl <- tab()
      #if these 2 values are changed, the filters are re-rendered because they are reactive values

      if (n > 0) {
        if(input$tableType=="Gene Expression"){
          isolate({lapply(seq_len(n), function(i) {
            #browser()
            list(
              selectInput(inputId = paste0("selctin", i),
                          label = paste0(i, ". Filter by:"),
                          choices = colnames(tbl[[2]])[1:2],
                          selected = AllInputs()[[paste0("selctin", i)]]), #call to AllInputs makes sure the user entered value stays constant through re-rendering
              textInput(inputId = paste0("textin", i),
                        label = paste0("Value:"),
                        value = AllInputs()[[paste0("textin", i)]])
            )
          })
          })
        }else{
          isolate({lapply(seq_len(n), function(i) {
            #browser()
            list( #putting these in a list makes sure they come out in the correct order: select1, text1, select2, text2, etc.
              selectInput(inputId = paste0("selctin", i),
                          label = paste0(i, ". Filter by:"),
                          choices = colnames(tbl[[2]]),
                          selected = AllInputs()[[paste0("selctin", i)]]),
              textInput(inputId = paste0("textin", i),
                          label = paste0("Value:"),
                          value = AllInputs()[[paste0("textin", i)]])
            )
          })
        })
        }
      }
    })

    #output UI list
    output$textbox_ui1 <- renderUI({ textboxes1() })

    #create "page" input
    output$length <- renderUI({
        tbl <- tab()
        #browser()
        selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    #create data table
    output$table <- DT::renderDataTable({

      tbl <- tab()
      num <- counter2$n

      #browser()

      if(num > 0){
        #searches for selct with a number on the end and gets all the inputs from inputId's like this
        y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
        y <- match(y,colnames(tbl[[2]]))
        y <- y[!is.na(y)]
        y <- data.frame(tbl[[2]][,y])

        z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
        #browser()
        x <- which(z!="",arr.ind = TRUE) #ignores all empty inputs for filters
        z <- z[x]
        #browser()

        if(identical(z,character(0))){
          z <- 0
        }else if(length(x) > 0){ #this statement makes sure the table is only being filtered by existing inputs
          z <- data.frame(strsplit(z,", ",fixed = TRUE)[[1]])
          z <- row.match(z, data.frame(y[,x])) #could make this more generalized and complicated with grep
        }else{
          z <- data.frame(strsplit(z,", ",fixed = TRUE)[[1]])
          z <- row.match(z, y) #could make this more generalized and complicated with grep
        }
      }else{
        z <- 0
      }

      #browser()

      #code to change "page" or dataset
      if(!is.data.frame(tbl[[2]])){
          n <- input$page
          n <- match(n, tbl[[1]]) #tbl[[1]] is a vector of data set names passed by the extraction functions
          #browser()
          if(is.null(n)){
            n <- 1
          }
          #browser()
          tbl <- tbl[[2]][[n]]
      }else{
        tbl <- tbl[[2]]
      }

      #this narrows down the table if the filters exist
      if(num==0 | identical(z,0) | identical(z,integer(0)))  {
        tbl
      }else{
        tbl[z,]
      }
    }, rownames = FALSE)
  })

  observe({
    input$DataID
    counter$countervalue <- 1 #sets to one whenever DataID is changed so that the user is asked which directory they want to save their data in later in the code, whenever this happens
  })

  counter <- reactiveValues(countervalue = 0,choice = getwd())

  #SAVE function
  observeEvent(input$save, {

    #browser()

    #set correct directory
    if(counter$countervalue == 1){
      if (exists('utils::choose.dir')) {
        counter$choice <- choose.dir(caption = "Select folder") #this is only available on Windows
      } else {
        counter$choice <- tk_choose.dir(caption = "Select folder")
      }
      setwd(counter$choice)
        #only call choose.dir if input$DataID is changed after save is hit
      d <- paste("data.",input$page,sep = "") #sets name of the new folder to data. the "page" or dataset selected
      dir.create(d)
      setwd(d)
    }else if(counter$countervalue > 1){ #if the data set hasn't been changed, the directory stays the same
      setwd(counter$choice)
      d <- paste("data.",input$page,sep = "")
      dir.create(d)
      setwd(d)
    }

    counter$countervalue <- counter$countervalue + 1

    #browser()

    #similar code to above, but has to keep the inputs from the text inputs, so that they can be used in the save names
    tbl <- tab()
    y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
    y <- match(y,colnames(tbl[[2]]))
    y <- y[!is.na(y)]
    y <- data.frame(tbl[[2]][,y])

    #browser()

    if(dim(y)[2] > 0){
      z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
      #browser()
      x <- which(z!="",arr.ind = TRUE)
      z <- z[x]
      #browser()

      if(identical(z,character(0))){
        z <- 0
      }else if(length(x) > 0){
        z <- data.frame(strsplit(z,", ",fixed = TRUE)[[1]])
        x <- row.match(z, data.frame(y[,x])) #could make this more generalized and complicated with grep
      }else{
        z <- data.frame(strsplit(z,", ",fixed = TRUE)[[1]])
        x <- row.match(z, y) #could make this more generalized and complicated with grep
      }
    }else{
      z <- 0
    }

    #browser()

    if(identical(z,0) | identical(z,integer(0)))  {
      tbl <- tbl[[2]]
      if(input$tableType == "Gene Expression"){
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.long.rds", sep="")) #adds a .long when the long data checkbox has been marked
        }else{
          saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
      }else{
        #browser()
        saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
      }
      #the next save functions take the vector of user text inputs and add them to the end of the file name
    }else{
      tbl <- tbl[[2]][x,]
      if(input$tableType == "Gene Expression"){
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.long.",paste(z[,1],collapse = "."),".rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(z[,1],collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        saveRDS(tbl, file = paste("geneAnnotation.",paste(z[,1],collapse = "."),".rds", sep=""))
      }else{
        #browser()
        saveRDS(tbl, file = paste("sampleAnnotation.",paste(z[,1],collapse = "."),".rds", sep=""))
      }
    }
  })
}
shinyApp(ui = ui, server = server)
