'''
Created on 24 Mar 2011

@author: Simon Bull
'''

import os
import sys
import random

def main(inFASTA, alignedHuman, length, toGenerate, resultsDir, makeInput=False, entireAlign=False):
    """inFASTA is the fasta file from which to generate the new fasta files
    alignedHuman is the alignment file of the entire human proteome
    length is the length of the fasta files to generate, and toGenerate is the number of fasta files to generate (both are lists)
    resultsDir is the directory to store the created files in
    makeInput means just use the sequences from inFASTA

    toGenerate and length match up in the following way:
    length = [A, B, C]  toGenerate = [X, Y, Z]
    then make X files containing A sequences, Y files containing B sequences and Z files containing C sequences
    A, B, C, X, Y and Z should all be integers
    """

    if entireAlign:
        filesCreated = []
        outputFile = resultsDir + '\\Aligned.txt'
        filesCreated.append(resultsDir + '\\Sequences.fasta')
        writeTo = open(outputFile, 'w')
        writeFasta = open(resultsDir + '\\Sequences.fasta', 'w')
        readFrom = open(inFASTA, 'r')
        for line in readFrom:
            writeFasta.write(line)
        readFrom.close()
        readFrom = open(alignedHuman, 'r')
        for line in readFrom:
            writeTo.write(line)
        readFrom.close()
        writeTo.close()
        writeFasta.close()
        return filesCreated

    readFrom = open(inFASTA, 'r')
    # Extract the proteins from the file and record them in listFasta
    listFasta = {}
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

    # Record all alignments between human proteins
    readFrom = open(alignedHuman, 'r')
    alignments = {}
    for line in readFrom:
        chunks = line.split()
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

    
    filesCreated = []

    if makeInput:
        outputFile = resultsDir + '\\Aligned.txt'
        filesCreated.append(resultsDir + '\\Sequences.fasta')

        writeTo = open(outputFile, 'w')
        writeFasta = open(resultsDir + '\\Sequences.fasta', 'w')
        for protein in listFasta.keys():
#            tempProtein = protein
            tempProtein = 'sp|' + protein  # Only use this for the stuff from the database or for files where the proteins don't begin with sp
            if alignments.has_key(tempProtein):
                for hit in alignments[tempProtein].keys():
#                    tempHit = hit
                    tempHit = hit[3:]  # Only use this for the stuff from the database or for files where the proteins don't begin with sp
                    if tempHit in listFasta.keys():
                        writeTo.write(alignments[tempProtein][hit])
            writeFasta.write('>' + tempProtein + '\n' + listFasta[protein] + '\n')
        writeTo.close()
        writeFasta.close()

        return filesCreated

    for i in range(len(length)):
        toGen = toGenerate[i]
            
        lengthOfFastaFile = length[i]
        print 'Generating the files with ', lengthOfFastaFile, ' sequences.'

        for j in range(toGen):
            print 'Generating file number ', j, ' out of ', toGen - 1
            if j < 10:
                numberToGen = '0' + str(j)
            else:
                numberToGen = str(j)

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
            

    return filesCreated
