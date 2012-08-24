# Takes two arguments. The first is a string of comma delimited positive integers (e.g. '1,1,1,2,2,3').
# The first two numbers indicate an undirected edge, the 3rd and 4th numbers another edge etc.
# The second argument is the file location to write the saved image to.

# Import the igraph library and graphics
library(igraph)
require(tcltk)

# Define the global variable to indicate the graph has been finished and the button press procedure
saved <- FALSE
discarded <- FALSE
waitForUserSave <- function()
{
	saved <<- TRUE
}
waitForUserDiscard <- function()
{
	discarded <<- TRUE
}

# Load in the arguments passed from the commandline
args <- commandArgs(TRUE)

# Split the string of integers into a vector of integers
params <- data.frame(strsplit(args,'-'))
outFile <- as.character(params[2,])
intArgs <- as.integer(unlist(strsplit(as.character(params[1,]),',')))

# Create the graph
g <- graph(c(intArgs ), directed=FALSE)

# Provide user information and create the saved button
'Please do not close the graph manually as this will effect the result. Please press the save or discard button when you are finished modifying the graph.'
tt <- tktoplevel()
saved.but <- tkbutton(tt, text = 'Save', command = waitForUserSave)
discard.but <- tkbutton(tt, text = 'Discard', command = waitForUserDiscard)
tkgrid(saved.but)
tkgrid(discard.but)
tkfocus(tt)

# Plot the graph
id = tkplot(g, layout = layout.fruchterman.reingold)

# Poll to see if user is finished
while (!saved && !discarded)
	Sys.sleep(2)

# Get user created coordinates
coords = tkplot.getcoords(id)

# Close the graph and button
tkplot.close(id)
tkdestroy(tt)

# Output the coordinates to a local file
# If the graph has three nodes 1, 2 and 3 with coords (20, 20), (30,100) and (50,90)
# then the text file will look like:
# 20
# 30
# 50
# 20
# 100
# 90
# i.e. x1 \n x2 \n ... xn \n y1 \n y2 \n ... yn
if (saved)
    write.graph(g, outFile, 'edgelist')
    #write(coords, outFile, 1, FALSE, '\t')