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

import processPISCES
import adjlistcreation

# Import the culling methods.
import CMBlastCuller
import CMVSA
import CMFIS
import CMGLP
import CMLeaf
import CMNeighbourCull
import CMUSearch


def main(args):
    """
    resultsDir is the directory where the results are
    algsToUse are the algorithms to use (names split by '-')
    timeitIterQuick is the number of times to run the Leaf algorithm
    """

    #####
    ## This assumes that all the BLASTing is done and the Aligned files from the BLASTing are in the resultsDir in their respective folders
    ## e.g. a crash has occured after all the BLASTing but the results are not generated yet so you need the results from the BLAST results
    ## Each folder must be cleared out of every file except aligned and the fasta file of sequences
    #####

    resultsDir, algsToUse, timeitIterQuick = args
    
    acceptedAlgs = ['Leaf', 'FIS', 'NeighbourCull', 'GLP', 'VSA', 'BlastCuller', 'UCLUST']
    tempAlgs = algsToUse.split('-')
    for i in tempAlgs:
        if i not in acceptedAlgs:
            print 'Algorithm ', i, ' is not recognised as a valid algorithm.'
            sys.exit()
    algorithmsToUse = [x for x in acceptedAlgs if x in tempAlgs]
    assert(len([i for i in algorithmsToUse if i in acceptedAlgs]) == len(algorithmsToUse))

    timeitIterQuick = int(timeitIterQuick)

    # Determine if the results directory is in fact a directory
    if not os.path.isdir(resultsDir):
        print 'The results directory is not a directory.'
        sys.exit()
    
    # Get the location of the PISCES bin directory
    currentFilePath = os.path.abspath( __file__ )
    PISCESPath = currentFilePath[::-1].split('\\',3)[-1][::-1] + '\\PISCES\\bin'
    if 'UCLUST' in algorithmsToUse:
        USearchLoc = currentFilePath[::-1].split('\\',3)[-1][::-1]
        USearch = USearchLoc + '\\usearch4.1.93_win32.exe'
        
    workingDirectory = os.getcwd()

    # Initialise the folders by removing everything except the Aligned file and the fasta file of sequences. Also set up the
    # dictionary with the locations for saving the results for each percentage.
    resultsFiles = dict([(i, resultsDir + '/' + i + 'PercentResults.txt') for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    dataDirs = os.listdir(resultsDir)
    for i in dataDirs:
        currentDir = resultsDir + '/' + i
        if os.path.isdir(currentDir):
            subDataDirs = os.listdir(currentDir)
            toKeep = [i + '.fasta', 'Aligned.txt']
            for j in subDataDirs:
                if j not in toKeep:
                    os.remove(currentDir + '/' + j)
        else:
            # Remove all non-directories.
            os.remove(resultsDir + '/' + i)
    dataDirs = os.listdir(resultsDir)

    for i in resultsFiles.keys():
        resultsFile = resultsFiles[i]
        resultsTabDelim = open(resultsFile, 'w')
        algNames = ''.join([i + '\t\t' for i in algorithmsToUse])
        headings = ''.join(['Time\tRemoved\t' for i in algorithmsToUse])
        resultsTabDelim.write('Test Set\tMean Degree\tNumber of Nodes\tNumber of Components\tLargest Component\tMean Component\tMean Degree of Largest Component\tPISCES\t\t' + algNames + '\n')
        resultsTabDelim.write('\t\t\t\t\t\t\tTime\tRemoved\t' + headings)

    # Go through dataDirs again. This time running the algorithms.
    for i in sorted(dataDirs):
        currentDir = resultsDir + '\\' + i
        print 'Now working on dataset ', i, ':'
        for percentage in sorted(resultsFiles.keys()):
            print '\tNow working on ', percentage, ' percent similarity cutoff'
            resultsTabDelim = open(resultsFiles[percentage], 'a')
            resultsTabDelim.write('\n' + i + '\t')

            print '\t\tNow timing PISCES'

            adjList, proteinNames = adjlistcreation.main(currentDir + '\\Aligned.txt', float(percentage))
            os.chdir(PISCESPath)
            os.rename(currentDir + '\\Aligned.txt', PISCESPath + '\\pdbaa.align')
            startPISCES = time.clock()
            subprocess.call('perl ' + PISCESPath + '\\Extract_Culled_SEQ.pl ' + percentage + ' 0 0 ' + currentDir + '\\' + i + '.fasta')
            PISCESTimed = time.clock() - startPISCES
            os.rename(PISCESPath + '\\log_pc' + percentage + '.log', currentDir + '\\' + percentage + 'PISCESCull.txt')
            # Process PISCES results to determine the number of proteins to remove and write out the time taken and number removed
            removedPISCES = processPISCES.main(currentDir + '\\' + percentage + 'PISCESCull.txt')
            writePISCESRemoved = open(currentDir + '\\' + percentage + 'PISCESCull.txt', 'w')
            for j in removedPISCES:
                writePISCESRemoved.write(j + '\n')
                writePISCESRemoved.close()
            os.rename(PISCESPath + '\\pdbaa.align', currentDir + '\\Aligned.txt')
            os.chdir(workingDirectory)

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
            degOut = open(currentDir + '\\' + percentage + 'DegreeHistogram.txt', 'w')
            degOut.write('Degree\tNumber of Occurences\n')
            for d in sorted(degree.keys()):
                degOut.write(str(d) + '\t' + str(degree[d]) + '\n')
            degOut.close()

            numOfNodes = len(tempAdjList.keys())
            resultsTabDelim.write(str(numOfNodes) + '\t')

            # Calculate component size information
            numberComponents = 0
            largestComponent = 0
            meanComponent = 0
            largestCompDeg = 0
            for comp in adjList.connectedcomponents():
                componentSize = len(comp)
                numberComponents += 1
                deg = [x for j in comp for x in tempAdjList[j]]
                deg = len(deg)
                if componentSize > largestComponent:
                    largestComponent = componentSize
                    largestCompDeg = float(deg) / largestComponent
                meanComponent += componentSize
            if numberComponents == 0:
                meanComponent = 0
            else:
                meanComponent = float(meanComponent) / numberComponents
            resultsTabDelim.write(str(numberComponents) + '\t')
            resultsTabDelim.write(str(largestComponent) + '\t')
            resultsTabDelim.write(str(meanComponent) + '\t')
            resultsTabDelim.write(str(largestCompDeg) + '\t')

            # Write out the PISCES results
            resultsTabDelim.write(str(PISCESTimed) + '\t' + str(len(removedPISCES)))

            # For each algorithm determine what it is and run the appropriate code
            multipleOfLeafTime = 10
            for alg in algorithmsToUse:
                if alg == 'Leaf':
                    print '\t\tNow using algorithm Leaf. ', time.clock()
                    LeafTimed = 0.0
                    for j in range(timeitIterQuick):
                        removedLeaf, proteinsToKeep, removeNode, nodesToKeep, timeTaken = CMLeaf.main(adjList, proteinNames)
                        LeafTimed = LeafTimed + timeTaken
                    LeafTimed /= float(timeitIterQuick)
                    LeafOutput = open(currentDir + '\\' + percentage + 'LeafCull.txt', 'w')
                    for j in removedLeaf:
                        LeafOutput.write(j + '\n')
                    LeafOutput.close()
                    resultsTabDelim.write('\t' + str(LeafTimed) + '\t' + str(len(removedLeaf)))

                    print '\t\tFinished algorithm Leaf. ', time.clock()
                    
                    timeAllowed = max(120, LeafTimed*multipleOfLeafTime)

                elif alg == 'FIS':
                    print '\t\tNow using algorithm FIS. ', time.clock()
                    removedFIS, proteinsToKeep, removeNode, nodesToKeep, FISTimed = CMFIS.main(adjList, proteinNames, timeAllowed)
                    FISOutput = open(currentDir + '\\' + percentage + 'FISCull.txt', 'w')
                    for j in removedFIS:
                        FISOutput.write(j + '\n')
                    FISOutput.close()
                    resultsTabDelim.write('\t' + str(FISTimed) + '\t' + str(len(removedFIS)))

                    print '\t\tFinished algorithm FIS. ', time.clock()
                    
                elif alg == 'NeighbourCull':
                    print '\t\tNow using algorithm NeighbourCull. ', time.clock()
                    if False:#percentage == '10':
                        resultsTabDelim.write('\t*\t*')
                    else:
                        removedNC, proteinsToKeep, removeNode, nodesToKeep, NCTimed, outOfTime = CMNeighbourCull.main(adjList, proteinNames, timeAllowed)
                        NCOutput = open(currentDir + '\\' + percentage + 'NCCull.txt', 'w')
                        for j in removedNC:
                            NCOutput.write(j + '\n')
                        NCOutput.close()
                        if outOfTime:
                            resultsTabDelim.write('\t' + str(NCTimed) + '\t*')
                        else:
                            resultsTabDelim.write('\t' + str(NCTimed) + '\t' + str(len(removedNC))

                    print '\t\tFinished algorithm NeighbourCull. ', time.clock()

                elif alg == 'VSA':
                    print '\t\tNow using algorithm VSA. ', time.clock()
                    removedVSA, proteinsToKeep, removeNode, nodesToKeep, VSATimed = CMVSA.main(adjList, proteinNames, timeAllowed)
                    VSAOutput = open(currentDir + '\\' + percentage + 'VSACull.txt', 'w')
                    for j in removedVSA:
                        VSAOutput.write(j + '\n')
                    VSAOutput.close()
                    resultsTabDelim.write('\t' + str(VSATimed) + '\t' + str(len(removedVSA))

                    print '\t\tFinished algorithm VSA. ', time.clock()

                elif alg == 'BlastCuller':
                    print '\t\tNow using algorithm BlastCuller. ', time.clock()
                    removedBlastCuller, proteinsToKeep, removeNode, nodesToKeep, BlastCullerTimed = CMBlastCuller.main(adjList, proteinNames, timeAllowed)
                    BlastCullerOutput = open(currentDir + '\\' + percentage + 'BlastCullerCull.txt', 'w')
                    for j in removedBlastCuller:
                        BlastCullerOutput.write(j + '\n')
                    BlastCullerOutput.close()
                    resultsTabDelim.write('\t' + str(BlastCullerTimed) + '\t' + str(len(removedBlastCuller))

                    print '\t\tFinished algorithm BlastCuller. ', time.clock()
                    
                elif alg == 'GLP':
                    print '\t\tNow using algorithm GLP. ', time.clock()
                    removedGLP, proteinsToKeep, removeNode, nodesToKeep, GLPTimed = CMGLP.main(adjList, proteinNames, timeAllowed)
                    GLPOutput = open(currentDir + '\\' + percentage + 'GLPCull.txt', 'w')
                    for j in removedGLP:
                        GLPOutput.write(j + '\n')
                    GLPOutput.close()
                    resultsTabDelim.write('\t' + str(GLPTimed) + '\t' + str(len(removedGLP))

                    print '\t\tFinished algorithm GLP. ', time.clock()

                elif alg == 'UCLUST':
                    sortedFile = currentDir + '\\' + percentage + 'UCLUSTSorted.fasta'
                    clusterFile = currentDir + '\\' + percentage + 'UCLUSTClustered.fasta'
                    sortCommand = USearch + ' --sort ' + currentDir + '\\' + i + '.fasta' + ' --output ' + sortedFile
                    clusterCommand = USearch + ' --cluster ' + sortedFile + ' --seedsout ' + clusterFile + ' --id 0.' + str(percentage)

                    stmt = """subprocess.call(sortCommand)
subprocess.call(clusterCommand)"""
                    USearchTimer = timeit.Timer(stmt, setup="import subprocess; sortCommand=%r; clusterCommand=%r;" %(sortCommand, clusterCommand))
                    USearchTimed = USearchTimer.timeit(1)/1

                    # Determine number of clusters and therefore the proteins to keep and to remove
                    removedUSearch = CMUSearch.main(currentDir + '\\' + i + '.fasta', clusterFile)
                    USearchOutput = open(currentDir + '\\' + percentage + 'UCLUSTCull.txt', 'w')
                    for j in removedUSearch:
                        USearchOutput.write(j + '\n')
                    USearchOutput.close()
                    
                    resultsTabDelim.write('\t' + str(USearchTimed) + '\t' + str(len(removedUSearch)))
            
            resultsTabDelim.close()

            # Clean up the PISCES bin directory
            for filename in glob.glob(PISCESPath + '\\*.align'):
                os.remove(filename)
            for filename in glob.glob(PISCESPath + '\\*.tmp'):
                os.remove(filename)
            for filename in glob.glob(PISCESPath + '\\cullseq_*'):
                os.remove(filename)
                


if __name__ == '__main__':
    main(sys.argv[1:])
