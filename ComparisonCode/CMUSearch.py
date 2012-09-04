def main(inFasta, clusterFasta):
    """

    @param inFasta: The FASTA fromat file containing the proteins in the dataset that UCLUST was run on.
    @type inFasta : string (file location)
    @param clusterFasta: The file containing the proteins that represent the clusters found by UCLUST.
    @type clusterFasta : string (file location)
    return @type: list
    return @use : The proteins removed from the dataset to make it non-redundant.

    """

    # Get all the names of the proteins.
    allProts = []
    readIn = open(inFasta, 'r')
    for i in readIn:
        if i[0] == '>':
            allProts.append(i[1:-1])
    readIn.close()

    # Get the names of the proteins kept.
    keptProts = []
    readIn = open(clusterFasta, 'r')
    for i in readIn:
        if i[0] == '>':
            keptProts.append(i[1:-1])
    readIn.close()

    # Determine the proteins that were not kept.
    notKept = [i for i in allProts if i not in keptProts]

    return notKept
