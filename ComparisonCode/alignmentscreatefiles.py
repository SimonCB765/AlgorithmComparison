'''
Created on 24 Mar 2011

@author: Simon Bull
'''

import os
import random
import re
import sys

def main(inFASTA, alignedFile, length, toGenerate, resultsDir):
    """Creates a set of datasets of user specified sizes, along with the alignment files for the datasets.

    Assumes that the FASTA file dataset has already been BLASTed, and the alignments determined.

    The size and number of the datasets generated is specified as follows:
    if length = [A, B, C] and toGenerate = [X, Y, Z]
    then X datasets of A sequences are produced, Y datasets of B sequences and Z datasets of C sequences

    @param inFASTA: The FASTA file from which the datasets will be generated.
    @type inFASTA : string (file location)
    @param alignedFile: The file of alignments of the proteins in the inFASTA file. Used to generate the alignment files for the generated datasets.
    @type alignedFile : string (file location)
    @param length: The number of proteins that should be in the datasets to be generated.
    @type length : list of integers
    @param toGenerate: The number of datasets of each size to generate.
    @type toGenerate : list of integers
    @param resultsDir: The directory where the directories containing the generated datasets and their alignment files should be stored.
    @type resultsDir : string (file location)

    """

    # Extract the proteins from the inFASTA and record them in listFasta.
    readFrom = open(inFASTA, 'r')
    listFasta = {}  # A dictionary mapping protein accessions to their sequences.
    currentProt = ''
    seq = ''
    first = True
    for line in readFrom:
        if line[0] == '>':
            if not first:
                listFasta[currentProt] = seq
                currentProt = ''
                seq = ''
            first = False
            currentProt = line[1:-1]
        else:
            seq += line[:-1]
    listFasta[currentProt] = seq
    readFrom.close()

    # Record all alignments between proteins.
    readFrom = open(alignedFile, 'r')
    alignments = {}
    for line in readFrom:
        # Correct the errors in the PISCES alignment file (if there are any). PISCES sometimes writes out alignments like this:
        # sp|Q8WZ42     34 11669   370-12005  sp|A6NNJ1    464     4   371    89    15  4.7e-002   388
        # The error is in the 370-12005 part. To fix this (in order that PISCES culling works properly) replace the '-' with a ' '.
        # PISCES will work even if the replacement is not made. However, there will be redundancy left in the 'non-redundant' dataset.
        hyphenSearch = re.search('[0,1,2,3,4,5,6,7,8,9]-[0,1,2,3,4,5,6,7,8,9]', line)
        if hyphenSearch:
            hyphenMatch = hyphenSearch.group(0)
            line = line.replace(hyphenMatch, hyphenMatch[0] + ' ' + hyphenMatch[2])
        chunks = line.split()
        query = chunks[0]
        hit = chunks[5]
        perc = float(chunks[10])
        if query == hit:
            continue
        if alignments.has_key(query):
            if alignments[query].has_key(hit):
                if alignments[query][hit][0] < perc:
                    alignments[query][hit] = [perc, line]
            else:
                alignments[query][hit] = [perc, line]
        else:
            alignments[query] = {}
            alignments[query][hit] = [perc, line]
        if alignments.has_key(hit):
            if alignments[hit].has_key(query):
                if alignments[hit][query][0] < perc:
                    alignments[hit][query] = [perc, line]
            else:
                alignments[hit][query] = [perc, line]
        else:
            alignments[hit] = {}
            alignments[hit][query] = [perc, line]
    readFrom.close()
    # Change the alignmentDicts dictionary so that it holds only strings, rather than the list of percent and string pairs.
    for i in alignments.keys():
        for j in alignments[i].keys():
            alignments[i][j] = alignments[i][j][1]

    for i in range(len(length)):
        toGen = toGenerate[i]
        lengthOfFastaFile = length[i]
        resultsSubDir = resultsDir + '\\Results' + str(lengthOfFastaFile)
        print 'Generating the datasets with ', lengthOfFastaFile, ' sequences.'

        for j in range(toGen):
            print 'Generating file number ', j, ' out of ', toGen - 1
            if j < 10:
                # Add a 0 on the front to get all the datasets having a two digit number (e.g. 0 -> 00).
                # Needs modifying if more than 99 datasets of a given size are being generated.
                numberToGen = '0' + str(j)
            else:
                numberToGen = str(j)

            # Create the directory that will hold the current dataset being generated.
            outputFolder = resultsSubDir + '\\' + str(lengthOfFastaFile) + '-' + numberToGen
            if os.path.isdir(outputFolder):
                confirm = raw_input(outputFolder + ' exists on the system. If you continue it will be deleted and recreated. Press \'y\' or \'Y\' without the \'\' to continue: ')
                if confirm.upper() == 'Y':
                    os.rmdir(outputFolder)
                    os.mkdir(outputFolder)
                else:
                    print 'Goodbye.'
                    sys.exit()
            elif not os.path.exists(outputFolder):
                os.mkdir(outputFolder)
            else:
                print outputFolder, ' exists and is not a directory. Please remove the file before retrying.'
                sys.exit()

            outputFile = outputFolder + '\\Aligned.txt'
            filesCreated.append(outputFolder + '\\' + str(lengthOfFastaFile) + '-' + numberToGen + '.fasta')

            # Choose a random subset of the entire set of proteins.
            chosen = random.sample(listFasta.keys(), lengthOfFastaFile)

            writeTo = open(outputFile, 'w')
            writeFasta = open(outputFolder + '\\' + str(lengthOfFastaFile) + '-' + numberToGen + '.fasta', 'w')
            for protein in chosen:
                if alignments.has_key(protein):
                    for hit in alignments[protein].keys():
                        if hit in chosen:
                            writeTo.write(alignments[protein][hit])
                writeFasta.write('>' + protein + '\n' + listFasta[protein] + '\n')
            writeTo.close()
            writeFasta.close()

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3].split('-'), sys.argv[4].split('-'), sys.argv[5])
