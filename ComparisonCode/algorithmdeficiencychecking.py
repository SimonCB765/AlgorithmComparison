import os
import shutil
import sys

def main(alignmentFile, resultsDir, algsToUse, outputDir):
    """Determines the level of known redundancy that algorithms failed to remove.

    Redundancy at a given percentage cutoff is deemed to be the number of proteins with a redundant relationship remaining in the non-redudnant dataset
    (i.e. the number of proteins that still have an incident edge remaining in the graph) divided by the number of proteins in the non-redundant dataset.

    Once culling results have been generated from a dataset, this script can be run to determine the level of redundancy (as
    calculated by PSI-BLAST) that remains in the non-redundant datasets.

    @param alignmentFile: The file containing all the possible alignments for the dataset that the sub-datasets (the datasets UCLUST was run on) were generated from.
    @type alignmentFile : string (file location)
    @param resultsDir: The directory containing the results of the culling. For more information on the structure of the directory see the README.
    @ytpe resultsDir : string (directory location)
    @param algsToUse: A list of the algorithms to include in the data generation, names separated by '-' (e.g. Leaf-GLP-VSA).
    @type algsToUse : string
    @param outputDir: The directory where the information about the redundancy which remains will be written.
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

    if os.path.isdir(outputDir):
        shutil.rmtree(outputDir)
    os.mkdir(outputDir)
    for i in percentages:
        writeTo = open(outputDir + '\\' + str(i) + '.txt', 'w')
        algNames = ''.join([i + '\t\t' for i in algsToUse])
        headings = ''.join(['\tProtsWithRedundantRelationship\tNonRedundantSetSize' for i in algsToUse])
        writeTo.write('Test Set\t' + algNames + '\n')
        writeTo.write(headings + '\n')
        writeTo.close()

    # Record all alignments between human proteins where the hit and query are not the same.
    readFrom = open(alignmentFile, 'r')
    for line in readFrom:
        chunks = line.split()
        if len(chunks) == 12:
            query = chunks[0]
            hit = chunks[4]
            perc = float(chunks[9])
        elif len(chunks) == 13:
            query = chunks[0]
            hit = chunks[5]
            perc = float(chunks[10])
        if query == hit:
            continue
        for i in percentages:
            if perc <= i:
                continue
            if alignmentDicts[i].has_key(query):
                if alignmentDicts[i][query].has_key(hit):
                    if alignmentDicts[i][query][hit] < perc:
                        alignmentDicts[i][query][hit] = perc
                else:
                    alignmentDicts[i][query][hit] = perc
            else:
                alignmentDicts[i][query] = {}
                alignmentDicts[i][query][hit] = perc
            if alignmentDicts[i].has_key(hit):
                if alignmentDicts[i][hit].has_key(query):
                    if alignmentDicts[i][hit][query] < perc:
                        alignmentDicts[i][hit][query] = perc
                else:
                    alignmentDicts[i][hit][query] = perc
            else:
                alignmentDicts[i][hit] = {}
                alignmentDicts[i][hit][query] = perc
    readFrom.close()

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
            algResults = []
            for k in algsToUse:
                # For each algorithm used, determine the proteins in the non-redundant set generated.
                errors = 0
                removedProteins = set([])
                currentFile = resultsDir + '\\' + i + '\\' + str(j) + k + 'Cull.txt'
                readIn = open(currentFile, 'r')
                for line in readIn:
                    removedProteins.add(line.strip())
                readIn.close()
                nonredundantSet = redundantDataset - removedProteins

                for p in nonredundantSet:
                    if not alignmentDicts[j].has_key(p):
                        continue
                    else:
                        similar = alignmentDicts[j][p].keys()
                        if len([h for h in similar if h in nonredundantSet]) > 0:
                            errors += 1
                algResults.append(str(errors))
                algResults.append(str(len(nonredundantSet)))

            # Write out the results.
            writeTo = open(outputDir + '\\' + str(j) + '.txt', 'a')
            writeTo.write(i + '\t' + '\t'.join(algResults) + '\n')
            writeTo.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
