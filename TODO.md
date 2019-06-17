[x] create base working package based on a main
[x] use Biobase to smooth out data table formation (going to break a lot of stuff in analyze used for graphs)
[x] Tidy up functions
[ ] refactor functions to be more general
    [x] ask if all geo data is in same format to implement already written functions
    [x] split up analysis and visualize into separate functions to prep data and graph/tabularize in different ways
    [x] update extraction functions to take a list of BioConductor formatted data, combine them, and look for a specific data set or extract a variable from the combined data
    [ ] update visual functions to take inputs and call the extraction functions to visualize specified sets of data
[ ]create shiny interface to search, tabularize, and graph different sets of data
    [ ] change the amount of inputs based on the type of data table selected
    [ ] create different titles for different types of data tables
    [ ] implement pre-defined graphs
    [ ] clean up column names for extExp to make it more readable in interface
