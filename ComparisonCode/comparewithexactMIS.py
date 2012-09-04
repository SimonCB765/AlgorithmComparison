import os
import shutil

def compare_exact_and_others(exactResults, otherAlgorithmsResults, outputLocation):
    otherAlgsUsed = otherAlgorithmsResults.keys()
    heading = 'Dataset\tExact Kept\tExact Time'
    for i in otherAlgsUsed:
        heading += ('\t' + i + ' Kept\t' + i + ' Kept Difference With Exact\t' + i + ' Time')
    results = dict([(i, dict([(j, [heading]) for j in ['20', '25', '30', '40', '50', '60', '70', '80', '90']]))
                    for i in ['100', '250', '500', '1000', '2000', '5000']])
    for cutoff in exactResults.keys():
        for dataset in sorted(exactResults[cutoff]):
            datasetSize = dataset.split('-')[0]
            line = (dataset + '\t' + str(exactResults[cutoff][dataset]['maxIndep']) + '\t' +
                    str(exactResults[cutoff][dataset]['time']))
            for i in otherAlgsUsed:
                otherAlgResults = otherAlgorithmsResults[i]
                if otherAlgResults[cutoff][dataset]['maxIndep'] == '*':
                    line += ('\t' + str(otherAlgResults[cutoff][dataset]['maxIndep']) + '\t*\t' + str(otherAlgResults[cutoff][dataset]['time']))
                else:
                    line += ('\t' + str(otherAlgResults[cutoff][dataset]['maxIndep']) + '\t' +
                             str(exactResults[cutoff][dataset]['maxIndep'] - otherAlgResults[cutoff][dataset]['maxIndep']) + '\t' +
                             str(otherAlgResults[cutoff][dataset]['time']))
            results[datasetSize][cutoff].append(line)

    if os.path.exists(outputLocation):
        # If it exists
        if os.path.isdir(outputLocation):
            # and is a directory then remove it.
            shutil.rmtree(outputLocation)
        else:
            # and it is not a directory then raise an error
            print 'The output location provided exists and is not a directory.'
            return        
    os.mkdir(outputLocation)
    for i in results.keys():
        outputDirectory = outputLocation + '\\\\' + i
        os.mkdir(outputDirectory)
        for j in results[i].keys():
            toWriteOut = '\n'.join(results[i][j])
            writeResults = open(outputDirectory + '\\\\' + i + '-' + j + '.txt', 'w')
            writeResults.write(toWriteOut)
            writeResults.close()

def process(folderWithMyResults):
    resultsPISCES = dict([(i, {}) for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    resultsLeaf = dict([(i, {}) for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    resultsFIS = dict([(i, {}) for i in ['10', '20', '35', '30', '40', '50', '60', '70', '80', '90']])
    resultsNeighbourCull = dict([(i, {}) for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    resultsVSA = dict([(i, {}) for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    resultsBlastCuller = dict([(i, {}) for i in ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']])
    for i in ['Results100', 'Results250', 'Results500', 'Results1000', 'Results2000', 'Results5000']:
        datasetSize = int(i[7:])
        currentDir = folderWithMyResults + '\\' + i
        dirs = os.listdir(currentDir)
        for j in dirs:
            if not os.path.isfile(currentDir + '\\' + j):
                continue
            cutoff = j[:2]
            readResults = open(currentDir + '\\' + j, 'r')
            readResults.readline()
            readResults.readline()
            for line in readResults:
                chunks = line.split()
                dataset = chunks[0].split('.')[0]
                PISCESTime = float(chunks[7])
                PISCESKept = datasetSize - int(chunks[8])
                leafTime = float(chunks[9])
                leafKept = datasetSize - int(chunks[10])
                FISTime = float(chunks[11])
                FISKept = datasetSize - int(chunks[12])
                neighbourCullTime = float(chunks[13])
                if chunks[14] == '*':
                    neighbourCullKept = chunks[14]
                else:
                    neighbourCullKept = datasetSize - int(chunks[14])
                VSATime = float(chunks[17])
                VSAKept = datasetSize - int(chunks[18])
                blastCullerTime = float(chunks[19])
                blastCullerKept = datasetSize - int(chunks[20])
                resultsPISCES[cutoff][dataset] = {}
                resultsPISCES[cutoff][dataset]['maxIndep'] = PISCESKept
                resultsPISCES[cutoff][dataset]['time'] = PISCESTime
                resultsLeaf[cutoff][dataset] = {}
                resultsLeaf[cutoff][dataset]['maxIndep'] = leafKept
                resultsLeaf[cutoff][dataset]['time'] = leafTime
                resultsFIS[cutoff][dataset] = {}
                resultsFIS[cutoff][dataset]['maxIndep'] = FISKept
                resultsFIS[cutoff][dataset]['time'] = FISTime
                resultsNeighbourCull[cutoff][dataset] = {}
                resultsNeighbourCull[cutoff][dataset]['maxIndep'] = neighbourCullKept
                resultsNeighbourCull[cutoff][dataset]['time'] = neighbourCullTime
                resultsVSA[cutoff][dataset] = {}
                resultsVSA[cutoff][dataset]['maxIndep'] = VSAKept
                resultsVSA[cutoff][dataset]['time'] = VSATime
                resultsBlastCuller[cutoff][dataset] = {}
                resultsBlastCuller[cutoff][dataset]['maxIndep'] = blastCullerKept
                resultsBlastCuller[cutoff][dataset]['time'] = blastCullerTime
            readResults.close()

    return resultsPISCES, resultsLeaf, resultsFIS, resultsNeighbourCull, resultsVSA, resultsBlastCuller

def processcsv(fileToCompare):
    results = dict([(i, {}) for i in ['20', '25', '30', '40', '50', '60', '70', '80', '90']])
    readResults = open(fileToCompare, 'r')
    readResults.readline()
    for line in readResults:
        line.strip()
        chunks = line.split(',')
        dataset = chunks[0]
        cutoff = chunks[1]
        maxIndep = chunks[3]
        time = chunks[4]
        results[cutoff][dataset] = {}
        results[cutoff][dataset]['maxIndep'] = int(maxIndep)
        results[cutoff][dataset]['time'] = float(time)
    readResults.close()

    return results

resultsPISCES, resultsLeaf, resultsFIS, resultsNeighbourCull, resultsVSA, resultsBlastCuller = process('G:\PhD-Snapshots\RedundancyRemovalResults\ResultsForExactComparison')
resultsExact = processcsv('G:\PhD-Snapshots\RedundancyRemovalResults\ResultsForExactComparison\ExactMisData\\xmips.csv')
resultsGLP = processcsv('G:\PhD-Snapshots\RedundancyRemovalResults\ResultsForExactComparison\ExactMisData\\mips.csv')

compare_exact_and_others(resultsExact,
                         {'PISCES' : resultsPISCES, 'Leaf' : resultsLeaf, 'FIS' : resultsFIS, 'NeighbourCull' : resultsNeighbourCull,
                          'GLP' : resultsGLP, 'VSA' : resultsVSA, 'BlastCuller' : resultsBlastCuller},
                         'G:\PhD-Snapshots\RedundancyRemovalResults\ResultsForExactComparison\ComparisonResults')
