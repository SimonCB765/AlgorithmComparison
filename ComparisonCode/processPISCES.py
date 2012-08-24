'''
Created on 28 Mar 2011

@author: Simon Bull
'''

def main(cullLoc):
    """Determine the proteins removed by PISCES.

    @param cullLoc: The location of the file containing the removed proteins.
    @type cullLoc : string
    return @type: list
    return @use : The proteins that PISCES removed to make the non-redundant dataset.

    """

    PISCES = open(cullLoc, 'r')

    removed = []

    for line in PISCES:
        chunks = line.split()
        removed.append(chunks[2])

    PISCES.close()

    removed = list(set(removed))

    return removed