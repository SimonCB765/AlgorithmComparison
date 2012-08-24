'''
Created on 28 Mar 2011

@author: Simon Bull
'''

import time
import subprocess

def pruneGraph(adjList, IDs):
    """The method by which the Leaf algorithm determines which nodes to remove from the graph.
    
    @param adjList: An adjacency list representation of the protein similarity graph.
    @type adjList : dictionary
    @param IDs: A list of the numerical indices of the nodes in the graph.
    @type IDs : list
    
    """
    
    removeList = []
    nodesInGraph = adjList.keys()
    if nodesInGraph == []:
        return removeList

    # Determine number of neighbours for each node.
    neighbours = {}
    for i in adjList.keys():
        numNeighbours = len(adjList[i])
        if neighbours.has_key(numNeighbours):
            neighbours[numNeighbours] |= set([i])
        else:
            neighbours[numNeighbours] = set([i])
    # Fill in the blank keys.
    neighbourKeys = neighbours.keys()
    for i in set(range(max(neighbourKeys))) - set(neighbourKeys):
        if neighbours.has_key(i):
            continue
        else:
            neighbours[i] = set([])
        
    
    while True:

        # Determine the maximum number of neighbours.
        maxNeighbours = max(neighbours.keys())
        while neighbours[maxNeighbours] == set([]):
            del neighbours[maxNeighbours]
            maxNeighbours = max(neighbours.keys())
        
        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        nClique = 1
        while nClique <= maxNeighbours:
            # Get the nodes with nClique neighbours.
            if neighbours.has_key(nClique):
                nodesOfInterest = neighbours[nClique]
                # For every node of interest see if the neighbours of the node are all neighbours of each other (i.e. a clique).
                for i in nodesOfInterest:
                    neighboursOfInterest = set(adjList[i])
                    while len(neighboursOfInterest) > 1:
                        toCheck = neighboursOfInterest.pop()
                        if neighboursOfInterest <= set(adjList[toCheck]):
                            continue
                        else:
                            break
                    else:
                        toRemove = adjList[i]
                        neighbours[nClique] -= set([i])
                        neighbours[0] |= set([i])
                        for i in toRemove:
                            removeList.append(i)
                            # Update the list of neighbours for each node that toRemove is adjacent to.
                            for j in adjList[i]:
                                numNeighbours = len(adjList[j])
                                neighbours[numNeighbours] -= set([j])
                                neighbours[numNeighbours - 1] |= set([j])
                                adjList[j].remove(i)
                            # Update the adjacency list to reflect the removal of to remove.
                            numNeighbours = len(adjList[i])
                            neighbours[numNeighbours] -= set([i])
                            neighbours[0] |= set([i])
                            adjList[i] = []
                        nClique = 1
                        break
                else:
                    # No clique found.
                    nClique += 1
            else:
                # No nodes with the desired number of neighbours.
                nClique += 1

        ########################################
        # Perform the NeighbourCull operation. #
        ########################################
        # Re-calculate this, as it may have changed since it was last calculated.
        maxNeighbours = max(neighbours.keys())
        while neighbours[maxNeighbours] == set([]):
            del neighbours[maxNeighbours]
            maxNeighbours = max(neighbours.keys())

        # If there are no nodes with neighbours then exit.
        if maxNeighbours == 0:
            return removeList

        # Get the IDs of the nodes with the max number of neighbours.
        nodesWithMaxNeighbours = list(neighbours[maxNeighbours])
        # If there is more than one node with the maximum number of neighbours determine which node to remove.
        if len(nodesWithMaxNeighbours) != 1:
            # Determine the number of neighbours for each node.
            extendedNeighbourhood = [adjList[x] + [x] for x in nodesWithMaxNeighbours]
            extendedNeighbourhood = [set().union(*[adjList[i] for i in a]) for a in extendedNeighbourhood]
            # Determine the size of each extended neighbourhood, and which nodes have the min size.
            sizes = [len(x) for x in extendedNeighbourhood]
            minSize = min(sizes)
            toRemove = nodesWithMaxNeighbours[sizes.index(minSize)]
        else:
            toRemove = nodesWithMaxNeighbours[0]
      
        removeList.append(toRemove)
        # Update the list of neighbours for each node that toRemove is adjacent to.
        for i in adjList[toRemove]:
            numNeighbours = len(adjList[i])
            neighbours[numNeighbours] -= set([i])
            neighbours[numNeighbours - 1] |= set([i])
            adjList[i].remove(toRemove)
        # Update the adjacency list to reflect the removal of to remove.
        adjList[toRemove] = []
        neighbours[maxNeighbours] -= set([toRemove])
        neighbours[0] |= set([toRemove])

def main(adj, names):
    """Use my heuristic method to calculate the maximum independent set.

    adj is the adjacency list for the graph
    names is a list of the names of each row in the matrix with names[0] corresponding to adjacent[0]
    """

    startTime = time.clock()
    rem = pruneGraph(adj.adjList(), range(len(names)))
    proteinsToCull = [names[x] for x in rem]
    proteinsToKeep = [names[x] for x in range(len(names)) if x not in rem]

    endTime = time.clock()

    return proteinsToCull, proteinsToKeep, rem, [x for x in range(len(names)) if x not in rem], endTime-startTime







##    # Determine the connected components of adj
##    subgraphs = adj.connectedcomponents()
##
##    # Create an adjacency list for each of the connected components, and record it along with the node ID for each
##    # of the nodes in the component.
##    subgraphMatrices = []
##    for i in subgraphs:
##        subSet = sorted(i)
##        subMat = adj.takesquare(subSet)
##        subMat = subMat.adjList()
##        subgraphMatrices.append((subMat, subSet))
##
##    startTime = time.clock()
##
##    # Non parallel method of getting results
##    removeNode = []
##    nodesToKeep = []
##    for i in subgraphMatrices:
##        subSetNodes = i[1]
##        rem = pruneGraph(i[0], range(len(subSetNodes)))
##        extendRemove = [subSetNodes[x] for x in range(len(subSetNodes)) if x in rem]
##        removeNode.extend(extendRemove)
##        extendKeep = [subSetNodes[x] for x in range(len(subSetNodes)) if x not in rem]
##        nodesToKeep.extend(extendKeep)
##
##    proteinsToCull = [names[x] for x in removeNode]
##    proteinsToKeep = [names[x] for x in nodesToKeep]
##
##    endTime = time.clock()
##
##    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, endTime-startTime
