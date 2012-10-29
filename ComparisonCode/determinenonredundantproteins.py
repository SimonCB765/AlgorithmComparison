import os
import sys

def main(resultsDir, algsToUse, outputDir):
    """Determines the proteins kept by each algorithm at each percentage cutoff.

    @param resultsDir: The directory containing the results of the culling. For more information on the structure of the directory see the README.
    @ytpe resultsDir : string (directory location)
    @param algsToUse: A list of the algorithms to include in the data generation, names separated by '-' (e.g. Leaf-GLP-VSA).
    @type algsToUse : string
    @param outputDir: The directory where the information about the proteins kept will be written.
    @type outputDir : string (directory location)

    """

    # Make a dictionary with one entry for each percentage cutoff used. For each cutoff, this will record all alignments that have a sequence
    # identity greater than the cutoff key value.
    percentages = [10, 20, 25, 30, 40, 50, 60, 70, 80, 90]
    alignmentDicts = dict([(i, {}) for i in percentages])
    acceptedAlgs = ['PISCES', 'Leaf', 'FIS', 'NC', 'GLP', 'VSA', 'BlastCuller', 'UCLUST']
    tempAlgs = algsToUse.split('-')
    algsToUse = [x for x in acceptedAlgs if x in tempAlgs]
    if (len([i for i in algsToUse if i in acceptedAlgs]) != len(algsToUse)):
        for i in tempAlgs:
            if i not in acceptedAlgs:
                print 'Algorithm ', i, ' is not recognised as a valid algorithm.'
        sys.exit()

    for i in os.listdir(resultsDir):
        if not os.path.isdir(resultsDir + '\\' + i):
            continue
        print 'Currently working on directory ', i

        # Determine the ID of every protein in the redundant dataset.
        readIn = open(resultsDir + '\\' + i + '\\' + i + '.fasta', 'r')
        redundantDataset = set()
        for line in readIn:
            if line[0] == '>':
                proteinID = line.split()[0][1:]
                redundantDataset.add(proteinID)
        readIn.close()

        for j in percentages:
            # For each sequence identity cutoff used:
            for k in algsToUse:
                # For each algorithm used, determine the proteins in the non-redundant set generated.
                removedProteins = set([])
                currentFile = resultsDir + '\\' + i + '\\' + str(j) + k + 'Cull.txt'
                readIn = open(currentFile, 'r')
                for line in readIn:
                    removedProteins.add(line.strip())
                readIn.close()
                nonredundantSet = redundantDataset - removedProteins

                # Write out the results.
                writeTo = open(outputDir + '\\' + str(k) + '-' + i  + '-' + str(j) + '.txt', 'w')
                writeTo.write('\n'.join(nonredundantSet))
                writeTo.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
