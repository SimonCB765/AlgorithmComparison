'''
Created on 28 Mar 2011

@author: Simon Bull
'''

import time

def pruneByNeighbours(adjList, IDs, timeAllowed, startTime):
    """same as pruneGraph"""
    
    removeList = []
    
    while True:

        if time.clock() - startTime > timeAllowed:
            return removeList, True
        
        # Determine number of neighbours for each node
        neighbours = [len(adjList[k]) for k in adjList.keys()]
        
        # If there are no nodes with neighbours then exit
        maxNeighbours = max(neighbours)
        if maxNeighbours == 0:
            return removeList, False
        
        # Get the IDs of the nodes with the max number of neighbours
        nodesWithMaxNeighbours = [x for x in range(len(neighbours)) if neighbours[x] == maxNeighbours]
        # If there is more than one node with the maximum number of neighbours determine which node to remove
        if len(nodesWithMaxNeighbours) != 1:
            # Determine the number of neighbours for each node
            extendedNeighbourhood = [adjList[x] + [x] for x in nodesWithMaxNeighbours]
            extendedNeighbourhood = [set([x for i in a for x in adjList[i]]) for a in extendedNeighbourhood]
            # Determine the size of each extended neighbourhood, and which nodes have this max size
            sizes = [len(x) for x in extendedNeighbourhood]
            minSize = min(sizes)
            nodesWithMaxNeighbours = [nodesWithMaxNeighbours[x] for x in range(len(sizes)) if sizes[x] == minSize]
            
        toRemove = nodesWithMaxNeighbours[0]            
        removeList.append(toRemove)            
        # Update the list of neighbours for each node that toRemove is adjacent to
        for i in adjList[toRemove]:
            adjList[i].remove(toRemove)            
        # Update the adjacency list to reflect the removal of to remove
        adjList[toRemove] = []


def main(adj, names, timeAllowed):

    # Determine the connected components of adj
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

    # Non parallel method of getting results
    removeNode = []
    nodesToKeep = []
    outOfTime = False
    for i in subgraphMatrices:
        subSetNodes = i[1]
        rem, outOfTime = pruneByNeighbours(i[0], range(len(subSetNodes)), timeAllowed, startTime)
        extendRemove = [subSetNodes[x] for x in range(len(subSetNodes)) if x in rem]
        removeNode.extend(extendRemove)
        extendKeep = [subSetNodes[x] for x in range(len(subSetNodes)) if x not in rem]
        nodesToKeep.extend(extendKeep)
        if outOfTime:
            break

    proteinsToCull = [names[x] for x in removeNode]
    proteinsToKeep = [names[x] for x in nodesToKeep]

    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, time.clock()-startTime, outOfTime
