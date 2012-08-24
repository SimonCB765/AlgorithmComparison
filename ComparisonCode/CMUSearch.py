def main(inFasta, clusterFasta):

    # Get all the names of the proteins
    allProts = []
    readIn = open(inFasta, 'r')
    for i in readIn:
        if i[0] == '>':
            allProts.append(i[1:-1])
    readIn.close()

    # Get the names of the proteins kept
    keptProts = []
    readIn = open(clusterFasta, 'r')
    for i in readIn:
        if i[0] == '>':
            keptProts.append(i[1:-1])
    readIn.close()

    # Determine the proteins that were not kept
    notKept = [i for i in allProts if i not in keptProts]

    return notKept
