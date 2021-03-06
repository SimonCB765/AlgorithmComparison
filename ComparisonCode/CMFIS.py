'''
Created on 7 Feb 2011

@author: Simon Bull
'''

import time

def add_nodes(adjList, start, timeAllowed, startTime):
    """Create a maximal independent set from a starting subset of the nodes in the graph by adding new nodes.
    
    If we assume the graph under consideration to be G, then start is a subset of the nodes of G such that it forms
    an independent set. The goal is to add as many nodes as possible to the initial set start in order to build a
    maximal independent set that is as large as possible. This is done by repeatedly adding the node which will cause
    the fewest additional nodes to become adjacent to the independent set under consideration.
    
    @param adjList: An adjacency list representation of the graph
    @type adjList: dictionary of lists
    @param start: The initial set of nodes which is to be augmented by adding nodes to it.
    @type start: list
    @param timeAllowed: The number of seconds that the algorithm is allowed to run for.
    @type timeAllowed : float
    @param startTime: The time when the algorithm was started.
    @type startTime : float
    return @type: list, boolean
    return @use : The nodes form the maximal independent set found, whether the time limit was exceeded.
    
    """

    newStart = list(start)

    # Determine the nodes in the graph which are adjacent to the nodes in start
    adjNodes = [x for i in newStart for x in adjList[i]]
    adjNodes = set(adjNodes)
    # Determine the nodes which are not adjacent to the nodes in start
    nonAdjNodes = set(adjList.keys())

    while True:
        
        # Update the record of the nodes which are not adjacent to the nodes in adjNodes.
        nonAdjNodes.difference_update(adjNodes)
        
        # If nonAdjNodes is == [], then a maximal independent set has been found as no new nodes can be added to adjNodes
        if not nonAdjNodes:
            break

        # Find the node which if added to adjNodes will cause adjNodes to gain the fewest nodes.
        numAdded = 1000000  # Arbitrarily large number to start with.
        nodeToAdd = None
        newAdjNodes = None
        # Loop through the non-adjacent nodes and determine the number of nodes that each node in nonAdjNodes will cause
        # to be added to adjNodes.
        for i in nonAdjNodes:
            intersect = set(adjList[i]).intersection(nonAdjNodes)
            interLength = len(intersect)
            # If a new minimum has been found mark the node i as the new best node to add to adjNodes.
            if interLength < numAdded:
                if interLength == 1:
                    nodeToAdd = i
                    newAdjNodes = intersect
                    break
                numAdded = interLength
                nodeToAdd = i
                newAdjNodes = intersect
        # Add the node found to the recorded independent set.
        newStart.append(nodeToAdd)
        # Add to adjNodes all the new nodes that this node will cause to become adjacent to the independent set.
        adjNodes.update(newAdjNodes)

        if time.clock() - startTime > timeAllowed:
            return newStart, True

    return newStart, False


def swap_nodes(adjList, start, maxSetSize, timeAllowed, startTime):
    """Swap nodes that are not in the maximal independent set with those that are to try and increase the size.
    
    Works by finding a node i that is not in the maximal independent set that is adjacent to only one node j that is
    in the maximal independent set. i is then swapped for j, and the set is again made maximal by adding as many nodes
    as possible using the add_nodes function.
    
    @param adjList: The adjacency list representation of the graph in question
    @type adjList: dictionary of lists
    @param start: The maximal independent set to permute
    @type start: list
    @param maxSetSize: The size of the largest maximal independent set found so far
    @type maxSetSize: integer
    @param timeAllowed: The number of seconds that the algorithm is allowed to run for.
    @type timeAllowed : float
    @param startTime: The time when the algorithm was started.
    @type startTime : float
    return @type: list, integer, boolean
    return @use : The nodes form the maximal independent set found, the size of the maximal independent set found, whether the time limit was exceeded.
    
    """

    newStart = list(start)
    maxSet = start
    # If there is no change in the maximum size of the independent set found over one loop of all nodes then exit.
    changed = False

    while changed:
        changed = False

        # For each node in the adjacency list:
        for key in adjList:
            # Determine how many nodes it is adjacenct to that are in the maximal independent set.
            inter = list(set(adjList[key]).intersection(newStart))
            if len(inter) == 1:
                # If it is only adjacent to one then swap it with that one.
                temp = list(newStart)
                temp.remove(inter[0])
                temp.append(key)

                # Attempt to add nodes to the set
                newStart = add_nodes(adjList, temp)

                # If a new largest set has been found then record the set and its size
                if len(newStart) > maxSetSize:
                    maxSetSize = len(newStart)
                    maxSet = newStart
                    changed = True

            if time.clock() - startTime > timeAllowed:
                return maxSet, maxSetSize, True

    return maxSet, maxSetSize, False


def fish(adj, timeAllowed, startTime):
    """Finds a maximal independent set as an approximation to the maximum independent set.
    
    @param adj: An adjacency list which represents the graph in which the maximal independent set is to be found.
    @type adj: dictionary of lists
    @param timeAllowed: The number of seconds that the algorithm is allowed to run for.
    @type timeAllowed : float
    @param startTime: The time when the algorithm was started.
    @type startTime : float
    return @type: list, boolean
    return @use : The nodes form the maximal independent set found, whether the time limit was exceeded.
    
    """

    # Create an adjacency list where the entry in the dictionary adjList for node n will be list which contains as its
    # elements the nodes that make up the closed neighbourood of n, i.e. n and all the nodes it is connected to.
    adjList = {}
    for i in adj:
        inc = list(adj[i])
        inc.append(i)
        adjList[i] = inc

    # Determine which node in the graph has connections ot the fewest other nodes (this will be used as the start node
    # to grow the maximal independent set).
    numList = [len(adj[x]) for x in adj]
    leastConnections = min(numList)
    leastConnected = numList.index(leastConnections)

    maxSet = []
    maxSetSize = 0
    
    # Determine an initial maximal independent set by adding nodes to the set containing the least connected node.
    indSet, outOfTime = add_nodes(adjList, [leastConnected], timeAllowed, startTime)
    maxSet = list(indSet)
    maxSetSize = len(indSet)
    indSet.sort()

    if outOfTime:
        return maxSet, True

    # Alter the maximal independent set found by swapping nodes in the set for nodes not in it in an attempt
    # to increase the number of nodes to be kept.
    indSet, indSetSize, outOfTime = swap_nodes(adjList, indSet, maxSetSize, timeAllowed, startTime)
    if indSetSize > maxSetSize:
        maxSet = list(indSet)
        maxSetSize = indSetSize

    if outOfTime:
        return maxSet, True
    else:
        return maxSet, False

def main(adj, names, timeAllowed):
    """Use the FIS heuristic method to calculate an approximation to the maximum independent set.

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

    # Determine the IDs of the nodes to keep by running the fish algorithm, and from this determine the names of
    # the proteins to keep and remove.
    removeNode = []
    nodesToKeep = []
    for i in subgraphMatrices:
        subSetNodes = i[1]
        keep, outOfTime = fish(i[0], timeAllowed, startTime)
        extendRemove = [subSetNodes[x] for x in range(len(subSetNodes)) if x not in keep]
        removeNode.extend(extendRemove)
        extendKeep = [subSetNodes[x] for x in range(len(subSetNodes)) if x in keep]
        nodesToKeep.extend(extendKeep)
        if outOfTime:
            break
    proteinsToCull = [names[x] for x in removeNode]
    proteinsToKeep = [names[x] for x in nodesToKeep]

    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, time.clock()-startTime
