'''
Created on 24 Mar 2011

@author: Simon Bull
'''

import os
import glob
import shutil
import sys
import subprocess
import timeit
import time

import BHOSLIB

# Import the culling methods.
import CMBlastCuller
import CMBSK
import CMFIS
import CMGLP
import CMLeaf
import CMNeighbourCull




def main(args):
    """
    resultsDir is the directory containing the BHOSLIB benchmark descriptions
    algorithmsToUse is a list of the algorithms to use in the comparison
    timeitIter is the number of iterations for timing
    """

    resultsDir, algsToUse, timeitIterQuick = args
    
    acceptedAlgs = ['Leaf', 'FIS', 'NeighbourCull', 'GLP', 'VSA', 'BlastCuller']
    tempAlgs = algsToUse.split('-')
    for i in tempAlgs:
        if i not in acceptedAlgs:
            raise
    algorithmsToUse = [x for x in acceptedAlgs if x in tempAlgs]
    assert(len([i for i in algorithmsToUse if i in acceptedAlgs]) == len(algorithmsToUse))

    timeitIterQuick = int(timeitIterQuick)

    toRun = os.listdir(resultsDir)
    for i in toRun:
        # Remove all the non-BHOSLIB graph files (i.e. the results files).
        if i[-4:] != '.mis':
            os.remove(resultsDir + '\\' + i)
    toRun = os.listdir(resultsDir)
    resultsFile = resultsDir + '\\Results.txt'
    resultsTabDelim = open(resultsFile, 'w')
    algNames = ''.join([i + '\t\t\t' for i in algorithmsToUse])
    headings = ''.join(['Time\tKept\tDifference to MIS\t' for i in algorithmsToUse])
    resultsTabDelim.write('Test Set\tMean Degree\tNumber of Nodes\t' + algNames + '\n')
    resultsTabDelim.write('\t\t\t' + headings + '\n')
    
    # Run the algorithms on each test file
    for benchFile in toRun:
        MISNum = benchFile.split('-')
        MISNum = int(MISNum[0][3:])
        benchmark = resultsDir + '\\' + benchFile
        bench = benchFile[:-4]
        print 'Now working on file ', bench
        resultsTabDelim.write(bench + '\t')
            
        adjList = BHOSLIB.main(benchmark)
            
        # Calculate the degree of all the nodes in the graph
        degree = {}
        meanDegree = 0
        tempAdjList = adjList.adjList()
        for key in tempAdjList.keys():
            deg = len(tempAdjList[key])
            if deg in degree:
                degree[deg] += 1
            else:
                degree[deg] = 1
            meanDegree += deg
        # If there are no nodes to remove then don't bother calculating the degree
        if len(tempAdjList.keys()) != 0:
            meanDegree /= float(len(tempAdjList.keys()))
        
        resultsTabDelim.write(str(meanDegree) + '\t')
        degOut = open(resultsDir + '\\' + bench + '-Histogram.txt', 'w')
        degOut.write('Degree\tNumber of Occurences\n')
        for i in sorted(degree.keys()):
            degOut.write(str(i) + '\t' + str(degree[i]) + '\n')
        degOut.close()

        numOfNodes = len(tempAdjList.keys())
        resultsTabDelim.write(str(numOfNodes))
        proteinNames = range(1, numOfNodes + 1)
            
        # For each algorithm determine what it is and run the appropriate code
        multipleOfLeafTime = 10
        for alg in algorithmsToUse:
                
            if alg == 'Leaf':
                print '\t\tNow using algorithm Leaf. ', time.clock()
                LeafTimed = 0.0
                for i in range(timeitIterQuick):
                    removedLeaf, proteinsToKeep, removeNode, nodesToKeep, timeTaken = CMLeaf.main(adjList, range(1, numOfNodes + 1))
                    LeafTimed = LeafTimed + timeTaken
                LeafTimed /= float(timeitIterQuick)
                LeafOutput = open(resultsDir + '\\' + bench + '-Leaf.txt', 'w')
                for i in proteinsToKeep:
                    LeafOutput.write(str(int(i) + 1) + '\n')
                LeafOutput.close()
                resultsTabDelim.write('\t' + str(LeafTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm Leaf. ', time.clock()
                
                timeAllowed = LeafTimed * multipleOfLeafTime
            
            elif alg == 'FIS':
                print '\t\tNow using algorithm FIS. ', time.clock()
                removedFIS, proteinsToKeep, removeNode, nodesToKeep, FISTimed = CMFIS.main(adjList, range(1, numOfNodes + 1), timeAllowed)
                FISOutput = open(resultsDir + '\\' + bench + '-FIS.txt', 'w')
                for i in proteinsToKeep:
                    FISOutput.write(str(int(i) + 1) + '\n')
                FISOutput.close()
                numRemoved = numOfNodes - len(nodesToKeep)
                resultsTabDelim.write('\t' + str(FISTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm FIS. ', time.clock()
                
            elif alg == 'NeighbourCull':
                print '\t\tNow using algorithm NeighbourCull. ', time.clock()
                removedNC, proteinsToKeep, removeNode, nodesToKeep, NCTimed, outOfTime = CMNeighbourCull.main(adjList, range(1, numOfNodes + 1), timeAllowed)
                NCOutput = open(resultsDir + '\\' + bench + '-NeighbourCull.txt', 'w')
                for i in removedNC:
                    NCOutput.write(str(int(i) + 1) + '\n')
                NCOutput.close()
                numRemoved = numOfNodes - len(nodesToKeep)
                if outOfTime:
                    resultsTabDelim.write('\t' + str(NCTimed) + '\t*\t*')
                else:
                    resultsTabDelim.write('\t' + str(NCTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm NeighbourCull. ', time.clock()

            elif alg == 'VSA':
                print '\t\tNow using algorithm VSA. ', time.clock()
                removedBSK, proteinsToKeep, removeNode, nodesToKeep, BSKTimed = CMBSK.main(adjList, range(1, numOfNodes + 1), timeAllowed)
                BSKOutput = open(resultsDir + '\\' + bench + '-VSA.txt', 'w')
                for i in removedBSK:
                    BSKOutput.write(str(int(i) + 1) + '\n')
                BSKOutput.close()
                numRemoved = numOfNodes - len(nodesToKeep)
                resultsTabDelim.write('\t' + str(BSKTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm VSA. ', time.clock()

            elif alg == 'BlastCuller':
                print '\t\tNow using algorithm BlastCuller. ', time.clock()
                removedBlastCuller, proteinsToKeep, removeNode, nodesToKeep, BlastCullerTimed = CMBlastCuller.main(adjList, range(1, numOfNodes + 1), timeAllowed)
                BlastCullerOutput = open(resultsDir + '\\' + bench + '-BlastCuller.txt', 'w')
                for i in removedBlastCuller:
                    BlastCullerOutput.write(str(int(i) + 1) + '\n')
                BlastCullerOutput.close()
                numRemoved = numOfNodes - len(nodesToKeep)
                resultsTabDelim.write('\t' + str(BlastCullerTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm BlastCuller. ', time.clock()
                
            elif alg == 'GLP':
                print '\t\tNow using algorithm GLP. ', time.clock()
                removedGLP, proteinsToKeep, removeNode, nodesToKeep, GLPTimed = CMGLP.main(adjList, range(1, numOfNodes + 1), timeAllowed)
                GLPOutput = open(resultsDir + '\\' + bench + '-GLP.txt', 'w')
                for i in removedGLP:
                    GLPOutput.write(str(int(i) + 1) + '\n')
                GLPOutput.close()
                numRemoved = numOfNodes - len(nodesToKeep)
                resultsTabDelim.write('\t' + str(GLPTimed) + '\t' + str(len(nodesToKeep)) + '\t' + str(MISNum - len(nodesToKeep)))

                print '\t\tFinished algorithm GLP. ', time.clock()
        
        resultsTabDelim.write('\n')
    
    resultsTabDelim.close()


if __name__ == '__main__':
    main(sys.argv[1:])
