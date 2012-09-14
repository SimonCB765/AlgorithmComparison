import os
import sys

def processdir(alignments, UCLUSTDir, percentage):
    """

    @param alignments: All the alignmnets with a percentage sequence identity greater than the cutoff percentage used to generate the non-redundant dataset.
    @type alignments : dictionary
    @param UCLUSTDir: A directory containing all the non-redundant datasets generated at a given percentage cutoff.
    @type UCLUSTDir : string (directory location)
    @param percentage: The percentage sequence identity used to generate the non-redundant dataset.
    @type percentage : integer
    return @type: list, list
    return @use : The number of proteins with a redudant relationship in each 'non-redundant' dataset at the percentage cutoff, the number of proteins in each 'non-redundant' dataset at the percentage cutoff.

    """
    
    numErrors = []  # The number of proteins with a redundant relationship remaining in each 'non-redundant' dataset.
    proteinLengths = []  # The number of proteins in each 'non-redundant' dataset.

    # Get the files in the results directory.
    filesUCLUST = os.listdir(UCLUSTDir)

    for f in filesUCLUST:
        print 'Now doing file: ', f, ' at percentage: ', percentage
        errors = 0
        proteins = []
        readIn = open(UCLUSTDir + '\\' + f, 'r')
        for line in readIn:
            if line[0] == '>':
                proteins.append(line[1:-1])
        readIn.close()
        proteinLengths.append(len(proteins))

        for p in proteins:
            if not alignments.has_key(p):
                continue
            else:
                similar = [i for i in alignments[p].keys() if alignments[p][i][10] > percentage]
                if len([i for i in similar if i in proteins]) > 0:
                    errors += 1
        numErrors.append(errors)

    return numErrors, proteinLengths

def main(alignmentFile, UCLUSTDir, outFile):
    """Determines the level of known redundancy that UCLUST failed to remove.

    Redundancy at a given percentage cutoff is deemed to be the number of proteins with a redundant relationship remaining in the dataset(s)
    divided by the number of proteins in the dataset(s). This gives the mean percentage of proteins that will have a redundant relationship
    that wasn't removed (i.e. the mean number of proteins that can expect to still have an edge remaining in the graph).

    Once UCLUST results have been generated from the sub-datasets of a dataset, this script can be run to determine the level of redundancy (as
    calculated by PSI-BLAST) that remains in the datasets after culling using UCLUST.

    @param alignmentFile: The file containing all the possible alignments for the dataset that the sub-datasets (the datasets UCLUST was run on) were generated from.
    @type alignmentFile : string (file location)
    @param UCLUSTDir: The directory containing the information about the UCLUST clusters generated. The directory must contain a subdirectory for
                      each percentage cutoff that was used, and that you want to generate remaining redundancy information about. Within each
                      percentage cutoff directory the results of the clustering for all of the sub-datasets must be placed. See the README for more
                      information on directory structure.
    @type UCLUSTDir : string (directory location)
    @param outFile: The file where the information about the redudnacy which remains will be written.
    @type outFile : string (file location)

    """

    # Make a dictionary with one entry for each percentage cutoff used. For each cutoff, this will record all alignments that have a sequence
    # identity greater than the cutoff key value.
    percentages = [10, 20, 25, 30, 40, 50, 60, 70, 80, 90]
    dicts = dict([(i, {}) for i in percentages])

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
            if dicts[i].has_key(query):
                dicts[i][query][hit] = line
            else:
                dicts[i][query] = {}
                dicts[i][query][hit] = line
            if dicts[i].has_key(hit):
                dicts[i][hit][query] = line
            else:
                dicts[i][hit] = {}
                dicts[i][hit][query] = line
    readFrom.close()

    dirs = os.listdir(UCLUSTDir)

    errors = {}  # A record of the number of redundant relationships remaining at each percentage cutoff.
    numProts = {}  # A record of the number of proteins kept at each percentage cutoff.

    # For each redundant dataset, compare the UCLUST 'non-redundant' dataset with the sequence identities from the file of alignments.
    for i in sorted(dirs, reverse=True):
        print 'Now working on percentage: ', i
        errors[i], numProts[i] = processdir(dicts[int(i)], UCLUSTDir + '\\' + i, int(i))
        totalErrors = float(sum(errors[i]))
        totalProts = sum(numProts[i])
        percErrors = (totalErrors / totalProts) * 100
        writeOut = open(outFile, 'a')
        writeOut.write(str(i) + ' percent.\n-------------------------------------\n')
        writeOut.write(str(totalErrors) + '\t' + str(totalProts) + '\t' + str(percErrors))
        writeOut.write('\n-------------------------------------\n')
        writeOut.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
