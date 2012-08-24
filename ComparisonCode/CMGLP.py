'''
Created on 29 Mar 2011

@author: Simon Bull
'''

import random
import time

def rule2(adjList, K):
    # Define q = 4 as used in the paper
    q = 4
    Sq = []
    for k in adjList.keys():
        neighbourhood = adjList[k]
        deltaI = len(set(K).difference(neighbourhood))
        if deltaI >= q:
            Sq.append(k)
    
    # Select the node to use to perturb K
    if Sq == []:
        init = random.choice(adjList.keys())
    else:
        init = random.choice(Sq)
    
    # Perturb K
    tempK = set(K).intersection(adjList[init])
    tempK = tempK.union([init])
    
    return list(tempK)

def calcC(adjList, K, p):
    
    nodesSatisfyingCriteria = []
    
    for k in adjList.keys():
        if len(set(K).difference(adjList[k])) == p and k not in K:
            nodesSatisfyingCriteria.append(k)
    
    return nodesSatisfyingCriteria

def GLP(adjList, timeAllowed, startTime):
    
    KBest = []
    K = []
    KPrime = []
    
    while time.clock() - startTime < timeAllowed: 
        # Perturb K
        K = rule2(adjList, K)
        
        # Equivalent to U in the paper
        bannedNodes = []
        
        while time.clock() - startTime < timeAllowed:
            CZeroK = calcC(adjList, K, 0)
            COneK = calcC(adjList, K, 1)
            
            if set(CZeroK).difference(bannedNodes) != set([]):
                nodeI = random.choice(list(set(CZeroK).difference(bannedNodes)))
                K = list(set(K + [nodeI]))
                if bannedNodes == []:
                    KPrime = K
            elif set(COneK).difference(bannedNodes) != set([]):
                nodeI = random.choice(list(set(COneK).difference(bannedNodes)))
                setJ = set(K).difference(adjList[nodeI])
                K = set(K + [nodeI])
                K = list(K.difference(setJ))
                bannedNodes = list(set(bannedNodes).union(setJ))
            elif set(bannedNodes).intersection(CZeroK) != set([]):
                nodeI = random.choice(list(set(bannedNodes).intersection(CZeroK)))
                K = list(set(K + [nodeI]))
            
            CZeroK = calcC(adjList, K, 0)
            COneK = calcC(adjList, K, 1)

            if ((CZeroK != [] or set(COneK).difference(bannedNodes) != set([])) and
                (set(KPrime).intersection(K) != set([]))):
                continue
            else:
                break
        
        if len(K) > len(KBest):
            KBest = K

    return KBest
       
        

def main(adj, names, timeAllowed):
    """Use the GLP heuristic method to calculate an approximation to the maximum independent set.

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

    complemenetSubgraphMatrices = []
    totalNodes = 0
    for i in subgraphMatrices:
        oldSubsetNodes = i[1]
        numberOldNodes = len(oldSubsetNodes)
        
        # Calculate the complement of the adjacency list as GLP works by finding cliques. Also prevent self loops in the complement
        # Record any nodes that are connected to all others as these will be isolated nodes in the complement
        indAdjList = i[0]        
        cliqueAdjList = {}
        isolated = []
        subsetNodes = []
        for k in indAdjList.keys():
            if len(indAdjList[k]) == numberOldNodes - 1:
                isolated.append(k)
            else:
                cliqueAdjList[k] = [x for x in range(numberOldNodes) if x not in indAdjList[k]and x != k]
                totalNodes += 1
                subsetNodes.append(k)

        if subsetNodes == []:
            complemenetSubgraphMatrices.append((isolated[0],isolated))
        else:
            complemenetSubgraphMatrices.append((cliqueAdjList,subsetNodes))
    totalNodes = float(totalNodes)

    allStart = time.clock()

    removeNode = []
    nodesToKeep = []
    for i in complemenetSubgraphMatrices:

        if type(i[0]) == int:
            keep = [i[0]]
        else:
            startTime = time.clock()
            subgraphTime = timeAllowed * (len(i[0].keys()) / totalNodes)
            keep = GLP(i[0], subgraphTime, startTime)

        extendRemove = [i[1][x] for x in range(len(i[1])) if x not in keep]
        removeNode.extend(extendRemove)
        extendKeep = [i[1][x] for x in range(len(i[1])) if x in keep]
        nodesToKeep.extend(extendKeep)

    proteinsToCull = [names[x] for x in removeNode]
    proteinsToKeep = [names[x] for x in nodesToKeep]

    return proteinsToCull, proteinsToKeep, removeNode, nodesToKeep, time.clock()-allStart
