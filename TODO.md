[x] create base working package based on a main
[x] use Biobase to smooth out data table formation (going to break a lot of stuff in analyze used for graphs)
[x] Tidy up functions
[ ] refactor functions to be more general
    [x] ask if all geo data is in same format to implement already written functions
    [x] split up analysis and visualize into separate functions to prep data and graph/tabularize in different ways
    [x] update extraction functions to take a list of BioConductor formatted data, combine them, and look for a specific data set or extract a variable from the combined data
    [ ] implement generalized graphing (take filtered table, and if less than 10 rows, create different plots)
    [x] change any $Symbol variable to grep search of column names because some data sets use gene symbol instead of just symbol
[ ] create shiny interface to search, tabularize, and graph different sets of data
    [x] change the amount of inputs based on the type of data table selected
    [x] create different titles for different types of data tables
    [x] implement pre-defined graphs (change box to ggplot)
    [x] give user ability to save data 
    [x] find way to include row names
    [x] give user ability to choose pages of a data set with multiple data frames 
    [x] improved moving through multiple data sets
    [x] pass names of data sets back when returning a list
    [x] use selectInput to choose which data set and what type of data set to save
    [x] add search by sample ID
    [ ] allow user to choose which variable they want to search by instead of just geneSymbol, potentially allow addition of more filters
        [x] limit amount of filters based off number of columns
            [x] take away filters when type of data table is changed
        [x] pass choices to selectize (maybe just text input)
        [x] fix reset of selection
        [x] fix save function
        [x] remove debug numbering
        [x] only use inputs that exist for filtering
        [x] make it so that geneExpression tables only have 2 options for filters
        [x] implement numeric filters
        [ ] make an option for a filter disappear from the list for additional filters/make identical filter by's compatible
    [ ] look at DT table filters 
    [x] give checkbox option for long format of data which is easier for graphing
    [x] add ID column to extExp or move gene symbols to rownames
[x] fix extraction functions to work with lists with multiple data sets as input
    [x] change if statements and output because some inputs cant have all their data sets combined
    [x] clean up column names for extExp (necessary to fix extExp function)
    [x] fix symbol column
    [x] check if other functions need a similar update
    [x] test extGene and extSample for multiple data sets
[ ] implement recount2 data
[ ] DOCUMENTATION
[ ] fix issue where GEOquery, stringr aren't importing - make biobase and them required

