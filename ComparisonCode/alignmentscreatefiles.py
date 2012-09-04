'''
Created on 24 Mar 2011

@author: Simon Bull
'''

import os
import sys
import random

def main(inFASTA, alignedFile, length, toGenerate, resultsDir, makeInput=False, entireAlign=False):
    """Creates a set of datasets of user specified sizes, along with the alignment files for the datasets.

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
    @param resultsDir: The directory where the generated datasets and their alignment files should be stored.
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
        chunks = line.split()
        # It's necessary to check the number of splits on the line due to the PISCES alignment writing procedure writing alignments to file badly.
        if len(chunks) == 12:
            query = chunks[0]
            hit = chunks[4]
        elif len(chunks) == 13:
            query = chunks[0]
            hit = chunks[5]
        if query == hit:
            continue
        if alignments.has_key(query):
            alignments[query][hit] = line
        else:
            alignments[query] = {}
            alignments[query][hit] = line
        if alignments.has_key(hit):
            alignments[hit][query] = line
        else:
            alignments[hit] = {}
            alignments[hit][query] = line
    readFrom.close()

    for i in range(len(length)):
        toGen = toGenerate[i]
        lengthOfFastaFile = length[i]
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
            outputFolder = resultsDir + '\\' + str(lengthOfFastaFile) + '-' + numberToGen
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
