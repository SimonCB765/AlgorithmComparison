import sys

def main(args):
    """Generates the data about the lengths of the proteins in the entire dataset, and in the non-redundant ones.

    @param args: The command line arguments for the script.
    @type args : list

    The three elements of args are:
    args[0]: The directory containing the results of the culling. For more information on the structure of the directory see the README.
    args[1]: The name of the results directory subdirectory that contains the files with the proteins kept (and the FASTA file of all proteins).
    args[2]: A list of the algorithms to include in the data generation, names separated by '-' (e.g. Leaf-GLP-VSA).

    """

    resultsDir, dirName, algsToUse = args

    # The percentage sequence identity cutoffs to generate data for.
    cutOffsUsed = ['10', '20', '25', '30', '40', '50', '60', '70', '80', '90']

    # Get the results files.
    resultsFiles = [resultsDir + '/' + i + 'PercentResults.txt' for i in cutOffsUsed]

    # Get the FASTA format file that the non-redundant datasets were generated from.
    fastaFile = resultsDir + '\\' + dirName + '\\' + dirName + '.fasta'

    # Determine the length of every protein in the redundant dataset.
    proteinLengths = {}
    readIn = open(fastaFile, 'r')
    firstProtein = True
    for line in readIn:
        if line[0] == '>':
            if not firstProtein:
                proteinLengths[proteinID] = str(proteinLength)
            proteinID = line.split()[0][1:]
            proteinLength = 0
            firstProtein = False
        else:
            line = line.strip()
            proteinLength += len(line)
    proteinLengths[proteinID] = str(proteinLength)
    readIn.close()
    allProteins = proteinLengths.keys()

    redundantDatasetLengths = sorted([proteinLengths[i] for i in proteinLengths.keys()])

    algsToUse = algsToUse.split('-')
    lengths = {}
    for i in cutOffsUsed:
        # For each sequence identity cutoff used:
        lengths[i] = {}
        lengths[i]['Entire Dataset'] = redundantDatasetLengths
        for j in algsToUse:
            # For each algorithm used, determine the proteins in the non-redundant set generated, along with their lengths.
            nonredundantLengths = []
            removedProteins = set([])
            currentFile = resultsDir + '\\' + dirName + '\\' + i + j + 'Cull.txt'
            readIn = open(currentFile, 'r')
            for line in readIn:
                removedProteins.add(line.strip())
            for k in allProteins:
                if not k in removedProteins:
                    nonredundantLengths.append(proteinLengths[k])
            readIn.close()
            # Pad the record of the lengths of the kept proteins to make writing out the lengths easier.
            nonredundantLengths = sorted(nonredundantLengths)
            nonredundantLengths.extend(['' for k in xrange(len(redundantDatasetLengths) - len(nonredundantLengths))])
            lengths[i][j] = nonredundantLengths
        writeTo = open(resultsDir + '/' + i + 'HistogramData.txt', 'w')
        datasetsUsed = lengths[i].keys()
        writeTo.write('\t'.join(datasetsUsed) + '\n')
        for j in range(len(redundantDatasetLengths)):
            writeTo.write('\t'.join([lengths[i][k][j] for k in datasetsUsed]) + '\n')
        writeTo.close()


if __name__ == '__main__':
    main(sys.argv[1:])
