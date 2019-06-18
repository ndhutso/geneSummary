[x] create base working package based on a main
[x] use Biobase to smooth out data table formation (going to break a lot of stuff in analyze used for graphs)
[x] Tidy up functions
[ ] refactor functions to be more general
    [x] ask if all geo data is in same format to implement already written functions
    [x] split up analysis and visualize into separate functions to prep data and graph/tabularize in different ways
    [x] update extraction functions to take a list of BioConductor formatted data, combine them, and look for a specific data set or extract a variable from the combined data
    [ ] update visual functions to take inputs and call the extraction functions to visualize specified sets of data
    [x] change any $Symbol variable to grep search of column names because some data sets use gene symbol instead of just symbol
[ ] create shiny interface to search, tabularize, and graph different sets of data
    [x] change the amount of inputs based on the type of data table selected
    [x] create different titles for different types of data tables
    [ ] implement pre-defined graphs (change box to ggplot)
    [ ] give user ability to save data 
    [x] find way to include row names
    [ ] give user ability to choose pages of a data set with multiple data frames 
[ ] extExp does not work when a list with multiple data sets is input without a data name input
    [x] change if statements and output because some inputs cant have all their data sets combined
    [x] clean up column names for extExp (necessary to fix extExp function)
    [ ] fix symbol column
    [ ] check if other functions need a similar update
