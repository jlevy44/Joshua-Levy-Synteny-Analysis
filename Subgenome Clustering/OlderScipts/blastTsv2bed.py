#!/usr/bin/env python

blastFiles = [] 
fileList = open('BLASTtsv.list.txt', 'r')
for file in fileList:
    file.strip()
    blastFiles.append(file)

def call2properBedgraph(blastFiles):
    """Takes a list of genotype files with only one column for pos and converts them to proper bedgraph format to be sorted"""
    for blastFile in blastFiles:
        f = blastFile.rstrip()
        print f
        outFileName = f[:f.rfind('.')]+'.bed3'
        inputFile = open(f, 'r')
        open(outFileName, 'w').close()
        outputFile = open(outFileName, 'w')
        for line in inputFile:
            lineInList = line.split()
            lineOutputList = [lineInList[1], int(lineInList[8]), int(lineInList[8])+1]
            outputFile.write('%s\t%d\t%d\n' % tuple(lineOutputList))
        inputFile.close()
        outputFile.close()

call2properBedgraph(blastFiles)