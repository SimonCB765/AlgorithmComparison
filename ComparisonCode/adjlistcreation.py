'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import sparsematrix

def main(similarities, cutoffPercent = 20, maxEValue = 1, minAlignLength = 20):
    """Create a sparse matrix from the processed PSI-BLAST output.
    
    @param similarities: The location of the processed PSI-BLAST output.
    @type similarities: string
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent:  integer
    @param maxEValue: The maximum permissible value for the E value of an alignment.
    @type maxEValue: float
    @param minAlignLength: The number of amino acids aligned in the query and the hit sequence
                           must be >= this value for the percentage similarity to be deemed significant.
    @type minAlignLength: integer
    
    """
    
    proteinNames = []  # Store the names of all the proteins found to be too similar to another protein
    similarProteins = []  # Store the pairs that are too similar
    
    readSimilarities = open(similarities, 'r')
    
    for line in readSimilarities:
        
        chunks = line.split()
        if len(chunks) == 12:
            query = chunks[0]
            hit = chunks[4]
            percentage = chunks[9]
        elif len(chunks) == 13:
            query = chunks[0]
            hit = chunks[5]
            percentage = chunks[10]
        
        # Ignore similarities where the query and the hit are the same, the percentage similarity is <= cutoffPercent,
        # the E value is > maxEValue and the length of the alignment is < minAlignLength.
        invalid = (query == hit or
            float(percentage) <= cutoffPercent
            )
        # If the similarity is valid record the proteins as being too similar.
        if not invalid:
            proteinNames.append(query)
            proteinNames.append(hit)
            similarProteins.append(tuple(sorted((query, hit))))
    
    readSimilarities.close()
    
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
    
    return adjacent, proteinNames
