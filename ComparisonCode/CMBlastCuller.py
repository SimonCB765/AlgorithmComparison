'''
Created on 17 May 2011

@author: Simon Bull
'''

import time

# http://ieeexplore.ieee.org/Xplore/login.jsp?url=http%3A%2F%2Fieeexplore.ieee.org%2Fiel5%2F5162127%2F5162128%2F05162176.pdf%3Ftp%3D%26arnumber%3D5162176%26punumber%3D5162127&authDecision=-203
# A Graph Theoretic Algorithm for Removing Redundant Protein Sequences 
# Pengfei Liu  Zhenbing Zeng  Ziliang Qian  Kaiyan Feng  Yudong Cai

def BlastCuller(adjList, timeAllowed, startTime):

    remove = []
    notfinished = True

    while notfinished:

        if time.clock() - startTime > timeAllowed:
            return remove, True

        # Calculate degree of all node
        degree = dict([(i, len(adjList[i])) for i in adjList.keys()])

        # Find node with largest degree
        nodeWithMaxDeg = None
        maxDeg = 0
        for i in degree.keys():
            if degree[i] > maxDeg:
                nodeWithMaxDeg = i
                maxDeg = degree[i]

        remove.append(nodeWithMaxDeg)

        # Remvoe all links to nodeWithMaxDeg
        removedNodeNeighbours = adjList[nodeWithMaxDeg]
        for i in removedNodeNeighbours:
            adjList[i].remove(nodeWithMaxDeg)
            adjList[nodeWithMaxDeg] = []

        # Check if there are any nodes with neighbours
        notfinished = False
        for i in adjList.keys():
            if adjList[i] != []:
                notfinished = True
                break

    return remove, False
    

def main(adj, names, timeAllowed):
    """Use the BlastCuller heuristic method to calculate an approximation to the maximum independent set.

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

    """
    
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

    # Determine the IDs of the nodes to keep by running the fish algorithm, and from this determine the names of
    # the nodes to keep and remove.
    removeNode = []
    nodesToKeep = []
    for i in subgraphMatrices:
        subSetNodes = i[1]
        remove, outOfTime = BlastCuller(i[0], timeAllowed, startTime)
        extendRemove = [subSetNodes[x] for x in range(len(subSetNodes)) if x in remove]
        removeNode.extend(extendRemove)
        extendKeep = [subSetNodes[x] for x in range(len(subSetNodes)) if x not in remove]
        nodesToKeep.extend(extendKeep)
        if outOfTime:
            break

    proteinsToCull = [names[x] for x in removeNode]
    proteinsToKeep = [names[x] for x in nodesToKeep]

    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, time.clock()-startTime
