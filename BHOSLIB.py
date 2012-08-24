import sparsematrix

def main(DIMACSFile):
    """Take a file in DIMACS format and turn it into a sparse matrix format."""

    toNodes = []
    fromNodes = []
    numNodes = 0

    readIn = open(DIMACSFile, 'r')

    for line in readIn:
        if line[0] == 'c':
            continue
        elif line[0] == 'p':
            chunks = line.split()
            numNodes = int(chunks[2])
        elif line[0] == 'e':
            chunks = line.split()
            fromNodes.append(int(chunks[1]) - 1)
            toNodes.append(int(chunks[2]) - 1)

    readIn.close()

    adjacent = sparsematrix.sparse_matrix(numNodes)
    adjacent.addlist(fromNodes, toNodes)
    adjacent.addlist(toNodes, fromNodes)

    return adjacent
