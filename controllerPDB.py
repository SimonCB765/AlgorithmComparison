import gzip
import os
import sys

import sparsematrix
import CMLeaf

###################
import time
###################

def main(i, locForResults, dirPISCES, duplicates, resolutionData, fastaSequences, alignmentFile, minLength=20, maxLength=10000):
    """ Compare the culling algorithms desired with the downloadable lists from PISCES.
    i is the file being done
    """

    resultsFile = open(locForResults, 'a')

    chunks = i.split('_')
    numberPISCESProteinsKept = int(chunks[-1][6:-3])  # The number of proteins kept by PISCES
    seqIdenThreshold = int(chunks[1][2:])
    resolutionLimit = float(chunks[2][3:])
    RFactorLimit = float(chunks[3][1:])
    if len(chunks) == 7:
        inclNotXray = True
        inclCAOnly = False
    elif len(chunks) == 8:
        inclNotXray = True
        inclCAOnly = True
    else:
        inclNotXray = False
        inclCAOnly = False
    
    ##################
    startTime = time.clock()
    ##################
    
    resolutionInfo = process_resolution(resolutionData)
##    print 'Resolution information processed.'
    alignments = process_alignments(alignmentFile, seqIdenThreshold)
##    print 'Alignment information processed.'
    redundant = process_redundant(duplicates)
##    print 'Redundant proteins processed.'
    lengths = process_fasta(fastaSequences)
##    print 'Lengths processed.'

    ##################
##    print 'Time taken to process files: ', time.clock() - startTime
    startTime = time.clock()
    ##################
    print 'Now working on file: ', i
    
    # Go through resolution and select the identifiers of the proteins that meet the structural requirements
    validProteins = []
    for p in lengths.keys():
        proteinResolution = resolutionInfo[p]['Resolution']
        proteinRFactor = resolutionInfo[p]['RFactor']
        proteinFreeRFactor = resolutionInfo[p]['FreeRFactor']
        if not inclNotXray:
            # If non-Xray structures are not being included
            if resolutionInfo[p]['Experiment'] != 'XRAY' or proteinResolution == 'NA':
                # If the protein structure was not determined by Xray or the resolution is not known
                continue
        else:
            if resolutionInfo[p]['Experiment'] != 'XRAY':
                # If non-Xray structures are being included and the protein structure was not determined by Xray
                proteinResolution = 'NA'
                proteinRFactor = 'NA'
                proteinFreeRFactor = 'NA'
        if not inclCAOnly and resolutionInfo[p]['CAOnly'] == 'yes':
            # If the structure is only CA, and that type of structure is not being considered
            continue
        
        if proteinResolution != 'NA' and float(proteinResolution) > resolutionLimit:
            continue
        if lengths[p] < minLength or lengths[p] > maxLength:
            continue
        if proteinRFactor == 'NA':
            if not inclNotXray:
                continue
        elif float(proteinRFactor) > RFactorLimit:
            continue
        
        validProteins.append(p)
    
##        print 'Valid Proteins (i.e. pdblen): ', len(validProteins)

    # Free up memory used by the resolution and length dicts
    del resolutionInfo
    del lengths
    
    # Go through the list of proteins that meet the structural requirements and remove all that can be
    # replaced by a representative sequence. Add the representative to the list if not already in the list.
    nonRedundantProteins = set([])
    for p in validProteins:
        if redundant.has_key(p):
            nonRedundantProteins.add(redundant[p])
        else:
            nonRedundantProteins.add(p)

    # Free up memory used by the redundant dict
    del redundant
    
    # Go through the alignments and record all alignments for the proteins that remain
    proteinNames = []
    similarProteins = []
    singletonProteins = []
    for p in nonRedundantProteins:
        singleton = True
        if alignments.has_key(p):
            # If there are alignments involving p then record the ones that involve other proteins in nonRedundantProteins
            for h in alignments[p]:
                if h in nonRedundantProteins and alignments[p][h] > seqIdenThreshold:
                    singleton = False
                    # If the two proteins are too similar then record this
                    proteinNames.append(p)
                    proteinNames.append(h)
                    similarProteins.append(tuple(sorted([p, h])))
        if singleton:
            singletonProteins.append(p)

    # Free up memory used by the alignments dict
    del alignments
    
    proteinNames = list(set(proteinNames))
    proteinNames.sort()
    similarProteins = list(set(similarProteins))
    indexDict = dict((proteinNames[x], x) for x in range(len(proteinNames)))

    # Create the sparse matrix
    adjacent = sparsematrix.sparse_matrix(len(proteinNames))
    xValues = [indexDict[x] for (x,y) in similarProteins]
    yValues = [indexDict[y] for (x,y) in similarProteins]
    adjacent.addlist(xValues, yValues)
    adjacent.addlist(yValues, xValues)

    ##############################################################################################################

    # Calculate the degree of all the nodes in the graph
    degree = {}
    meanDegree = 0
    tempAdjList = adjacent.adjList()
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
    maxDegree = max(degree.keys())
    numOfNodes = len(tempAdjList.keys())
    
    print '\tNumber of Nodes: ', numOfNodes
    print '\tMax Degree: ', maxDegree
    print '\tMean Degree: ', meanDegree
    ##############################################################################################################
    
    removedLeaf, proteinsToKeep, removeNode, nodesToKeep, timeTaken = CMLeaf.main(adjacent, proteinNames)
    
    numberLeafProteinsKept = len(nonRedundantProteins) - len(removedLeaf)
    percentageImprovementLeaf = (float((numberLeafProteinsKept - numberPISCESProteinsKept)) / numberPISCESProteinsKept) * 100
    
    ######################
    print '\tLeaf kept: ', numberLeafProteinsKept
    print '\tPISCES kept: ', numberPISCESProteinsKept
    print '\tPercentage improvement: ', percentageImprovementLeaf
    print '\tTime taken: ', time.clock() - startTime
    ######################

    resultsFile.write(i + '\t' + str(numberPISCESProteinsKept) + '\t' + str(numberLeafProteinsKept) + '\t' + str(percentageImprovementLeaf) + '\t' +
                      str(numOfNodes) + '\t' + str(maxDegree) + '\t' + str(meanDegree) + '\n')
    resultsFile.close()

def process_redundant(duplicates):
    redundant = {}
    
    readRedundant = open(duplicates, 'r')
    for line in readRedundant:
        chunks = line.split()
        if len(chunks) == 3 and chunks[1] == 'by':
            redundantProt = chunks[0]
            representativeProt = chunks[2]
            redundant[redundantProt] = representativeProt
    readRedundant.close()
    
    return redundant

def process_fasta(fastaSequences):
    lengths = {}
    
    readFasta = open(fastaSequences, 'r')
    for line in readFasta:
        if line[0] == '>':
            chunks = line.split()
            protein = chunks[0][1:]
            lengths[protein] = int(chunks[1])
    readFasta.close()
    
    return lengths

def process_resolution(resolutionData):
    resolutionInfo = {}
    
    readResolution = open(resolutionData, 'r')
    for line in readResolution:
        chunks = line.split()
        protein = chunks[0]
        resolutionInfo[protein] = {}
        resolutionInfo[protein]['Experiment'] = chunks[1]
        resolutionInfo[protein]['Resolution'] = chunks[2]
        resolutionInfo[protein]['RFactor'] = chunks[3]
        resolutionInfo[protein]['FreeRFactor'] = chunks[4]
        resolutionInfo[protein]['CAOnly'] = chunks[5]
    readResolution.close()
    
    return resolutionInfo

def process_alignments(alignmentFile, seqIdenThreshold):
    alignments = {}
    
    alignFile = open(alignmentFile, 'r')
    for line in alignFile:
        chunks = line.split()
        query = chunks[0]
        hit = chunks[1]
        if query == hit:
            continue
        identity = float(chunks[2])
        if identity <= seqIdenThreshold:
            continue
        
        if not alignments.has_key(query):
            alignments[query] = {}
            
        if not alignments.has_key(hit):
            alignments[hit] = {}
            alignments[hit][query] = 0
        elif not alignments[hit].has_key(query):
            alignments[hit][query] = 0
        
        # The similarity is set to be the maximum of the two possible similarities (i.e. [query][hit] and [hit][query]
        maxIdentity = max(identity, alignments[hit][query])
        alignments[query][hit] = maxIdentity
        alignments[hit][query] = maxIdentity
    alignFile.close()
    
    return alignments


if __name__ == '__main__':
    BLASTDBFolder = sys.argv[1]
    duplicates = BLASTDBFolder + '\MakeNonRedundant.log.pdb'
    resolutionData = BLASTDBFolder + '\\resolution.dat'
    fastaSequences = BLASTDBFolder + '\\pdbaa'
    alignmentFile = BLASTDBFolder + '\\pdbaa.align'

    
    culledListFolder = sys.argv[2]  # Expects this input to be a folder that contains all the gzipped culled PISCES lists.
    resultsFile = culledListFolder + '\Results.txt'

    alreadyDone = []
    if os.path.isfile(resultsFile):
        readResults = open(resultsFile, 'r')
        for line in readResults:
            chunks = line.split()
            if chunks[0] == 'File':
                continue
            else:
                alreadyDone.append(chunks[0])
        readResults.close()
    else:
        writeResults = open(resultsFile, 'w')
        writeResults.write('File\tPISCES kept\tLeaf kept\t% improvement\tNumber of Nodes\tMax Degree\tMean Degree\n')
        writeResults.close()
    
    culledListsPISCES = os.listdir(culledListFolder)
    stillToDo = [i for i in culledListsPISCES if 'gz' in i]
    print len(stillToDo)
    stillToDo = [i for i in stillToDo if i not in alreadyDone]
    print len(stillToDo)
    
    for i in stillToDo:
        main(i, resultsFile, culledListFolder, duplicates, resolutionData, fastaSequences, alignmentFile)
