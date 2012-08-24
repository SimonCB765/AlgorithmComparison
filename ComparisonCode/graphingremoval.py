import igraph
import subprocess
import os

def main(adjList, outputLoc, removed, nodeLabels=[], interactive=False):
    """Generate an image of the protein similarity graph that displays the non-redundant and removed proteins.

    In the final graph the white nodes are the ones removed while the black ones are the ones remaining

    @param adjList: The adjacency matrix used to represent the graph
    @type adjList : dictionary
    @param outputLoc: The location where the image of the graph should be saved.
    @type outputLoc : string
    @param removed: The nodes that are no longer in the graph. The ids of the nodes run from 0 to n.
    @type removed : List of integers identifying the removed nodes.
    @param nodeLabels: The labels for the nodes in the graph.
    @type nodeLabels : A list of strings. If the list is empty then no labels are used. If the list is not empty but the length is less than the number of nodes in the graph, then default labels are used.
    @param interactive: Whether the user wants to manually layout the nodes of the graph.
    @type interactive : boolean

    """


    # Define useful variables
    numberNodes = len(adjList)

    # Convert the adjacency matrix to the list of edges string format for R and the edge list needed for Python construction
    edges = ''
    edgeList = []
    for i in adjList.keys():
        connectedTo = adjList[i]
        for j in connectedTo:
            if i < j:
                edges += str(i) + ',' + str(j) + ','
                edgeList.append((i, j))
    # Remove the trailing ','
    edges = edges[:-1]

    # If R fine tuning is selected
    if interactive:
        # Call the R script and create the folder for the output of it
        cwd = os.getcwd()
        outputDir = cwd + '\Temp'
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        savedCoords = outputDir + '\Visualisation.txt'
        RLoc = 'C:\Program Files\R\R-2.12.1\\bin\Rscript.exe'
        RScriptLoc = cwd + '\\interactiveGraphing.R'
        subprocess.call(RLoc + ' ' + RScriptLoc + ' ' + edges + ' ' + savedCoords)

        # Get the file written by R
        readCoords = open(savedCoords, 'r')
        coords = []
        for i in readCoords:
            coords.append(float(i[:-1]))

        # Reconfigure the coords so that rather than being a long list of all x coords followed by all y coords they are in the format igraph wants
        igraphCoords = []
        count = 0
        coordPos = 0
        while count < numberNodes:
            igraphCoords.append([coords[coordPos]])
            coordPos += 1
            count += 1
        count = 0
        while count < numberNodes:
            igraphCoords[count].append(coords[coordPos])
            coordPos += 1
            count += 1

    # Create the graph
    g = igraph.Graph(edgeList)

    # Define which vertices are to be removed, and mark them on the graph
    vertices = range(len(adjList))
    markRemoved = []
    for v in vertices:
        if v in removed:
            markRemoved.append('r')
        else:
            markRemoved.append('k')
    g.vs['removed'] = markRemoved
    colorDict = {'r' : 'white', 'k' : 'black'}
    g.vs['color'] = [colorDict[remove] for remove in g.vs['removed']]

    # If R fine tuning is selected
    if interactive:
        # Define the node layout
        layout = igraph.Layout(igraphCoords)
    # If interactive mode is diabled use the default Fruchterman-Reingold force-directed algorithm to determine node position
    else:
        layout = g.layout('fr')

    # Determine what the node labels should be
    if nodeLabels == []:
        g.vs['label'] = [''] * numberNodes
    elif len(nodeLabels) < numberNodes:
        print 'The number of labels supplied is less than the number of nodes. Default labels will therefore be used.'
    else:
        g.vs['label'] = nodeLabels
        

    # Save the graph
    igraph.plot(g, target = outputLoc, layout = layout)
