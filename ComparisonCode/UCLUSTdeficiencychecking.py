import os
import sys

def processdir(alignments, UCLUSTDir, percentage):
    
    numErrors = []
    proteinLengths = []

    # Get the files in the UCLUST directory
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

    @param alignmentFile:
    @type alignmentFile :
    @param UCLUSTDir:
    @type UCLUSTDir :
    @param outFile:
    @type outFile : strin
    return @type:
    return @use :

    """

    percentages = [10, 20, 30, 40, 50, 60, 70, 80, 90]

    dicts = dict([(i, {}) for i in percentages])

    # Record all alignments between human proteins where the hit and query are not the same
    readFrom = open(alignmentFile, 'r')
    alignments = {}
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

    errors = {}
    numProts = {}

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

    return errors, numProts


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
