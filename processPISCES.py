'''
Created on 28 Mar 2011

@author: Simon Bull
'''

def main(cullLoc):
    """Processes the PISCES removal file.

    records and returns in a list format the proteins that PISCES has culled

    cullLoc should be the location of the removal file."""

    PISCES = open(cullLoc, 'r')

    removed = []

    for line in PISCES:
        chunks = line.split()
        removed.append(chunks[2])

    PISCES.close()

    removed = list(set(removed))

    return removed