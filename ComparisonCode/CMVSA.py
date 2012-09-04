'''
Created on 17 May 2011

@author: Simon Bull
'''

import time


#
#
#
# http://camo.ici.ro/journal/vol12/v12a9.pdf A Simple Algorithm to Optimize Maximum Independent Set
# S. Balaji, V. Swaminathanb, K. Kannanc
#
#
#

def BSK(adjList, timeAllowed, startTime):
    """

    @param adjList: The adjacency list of the current connected component being examined.
    @type adjList : dictionary
    @param timeAllowed: The number of seconds that the algorithm is allowed to run for.
    @type timeAllowed : float
    @param startTime: The time when the algorithm was started.
    @type startTime : float
    return @type: list, boolean
    return @use : The nodes that compose the minimum vertex cover, whether the time limit was exceeded.

    """

    MVC = []
    notfinished = True

    while notfinished:

        # Calculate degree of every node.
        degree = dict([(i, len(adjList[i])) for i in adjList.keys()])

        # Calculate the support of all nodes, and the largest support of all nodes.
        support = {}
        largestSupport = 0
        for i in adjList.keys():
            supportSize = 0
            for j in adjList[i]:
                supportSize += degree[j]
            if supportSize > largestSupport:
                largestSupport = supportSize
            support[i] = supportSize

        # Find all nodes with the largest support.
        nodesWithMaxSupport = []
        for i in support.keys():
            if support[i] == largestSupport:
                nodesWithMaxSupport.append(i)

        # Determine which node to add to teh vertex cover.
        if len(nodesWithMaxSupport) == 1:
            # If only one node has the largest support, then add it to the vertex cover.
            removedNode = nodesWithMaxSupport[0]
            MVC.append(nodesWithMaxSupport[0])
        else:
            # If multiple nodes have the largest vertex cover, then find the node with the maximum degree.
            maxDegree = 0
            nodeWithMaxDeg = []
            for i in nodesWithMaxSupport:
                nodeDeg = degree[i]
                if nodeDeg > maxDegree:
                    maxDegree = nodeDeg
                    nodeWithMaxDeg = i
            # Add the node wiht the maximum degree to the vertex cover.
            removedNode = nodeWithMaxDeg
            MVC.append(nodeWithMaxDeg)

        # Remvoe all edges between removedNode and other nodes.
        removedNodeNeighbours = adjList[removedNode]
        for i in removedNodeNeighbours:
            adjList[i].remove(removedNode)
            adjList[removedNode] = []

        # Check if there are any nodes with neighbours.
        notfinished = False
        for i in adjList.keys():
            if adjList[i] != []:
                # If there are still nodes with neighbours, then continue searching for the minimum vertex cover.
                notfinished = True
                break

        if time.clock() - startTime > timeAllowed:
            return MVC, True

    return MVC, False
    

def main(adj, names, timeAllowed):
    """Use the VSA heuristic method to calculate an approximation to the maximum independent set.

    Returns a list of the proteins to keep and a list of the proteins to cull. The list of proteins to keep only contains the
    names of the proteins in the protein similarity graph that should be kept. If there are any proteins that were not
    included in adj (for example proteins with no neighbours), then these will NOT be included in the list of proteins to keep.
    See the README for a more in depth description of this.
    
    @param adj: A sparsematrix representation of the protein similarity graph
    @type adj : sparsematrix
    @param names: A list of the names of the proteins in adj. Ordered such that the name of the protein represented by node i
                  in adj is located at names[i].
    @type names : list
    @param timeAllowed: The maximum number of seconds the algorithm is allowed to run for.
    @type timeAllowed : float
    return @type: list, list, list, list, float
    return @use : names of the proteins to cull, names of the proteins from the graph to keep, numerical IDs of the proteins to cull, numerical IDs of the protein from the graph to keep, time taken by the algorithm

    """
    
    # Determine the connected components of the graph.
    subgraphs = adj.connectedcomponents()

    # Create an adjacency list for each of the connected components, and record it along with the node ID for each
    # of the nodes in the component.
    subgraphMatrices = []
    for i in subgraphs:
        subSet = sorted(i)
        subMat = adj.takesquare(subSet)
        subMat = subMat.adjList()
        subgraphMatrices.append((subMat, subSet))

    startTime = time.clock()

    # Determine the IDs of the nodes to keep by running the VSA algorithm, and from this determine the names of
    # the proteins to keep and remove.
    removeNode = []
    nodesToKeep = []
    for i in subgraphMatrices:
        subSetNodes = i[1]
        remove, outOfTime = BSK(i[0], timeAllowed, startTime)
        extendRemove = [subSetNodes[x] for x in range(len(subSetNodes)) if x in remove]
        removeNode.extend(extendRemove)
        extendKeep = [subSetNodes[x] for x in range(len(subSetNodes)) if x not in remove]
        nodesToKeep.extend(extendKeep)
    proteinsToCull = [names[x] for x in removeNode]
    proteinsToKeep = [names[x] for x in nodesToKeep]

    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, time.clock()-startTime
