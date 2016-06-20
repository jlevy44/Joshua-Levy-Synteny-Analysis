#Chr1	phytozomev11	gene	2903	10817	.	+	.	ID=LOC_Os01g01010.MSUv7.0;Name=LOC_Os01g01010
#Chr1	phytozomev11	mRNA	2903	10817	.	+	.	ID=LOC_Os01g01010.1.MSUv7.0;Name=LOC_Os01g01010.1;pacid=33123311;longest=1;Parent=LOC_Os01g01010.MSUv7.0
#Pahal.C01314.1	Chr03	7436628	7436867	80	1	U	1007309

#gffFiles= ['t.PAC2_0.313.gff3',	't.PAC2_0.323.gff3', 't.PAC2_0.308.gff3',
           #'t.PAC2_0.314.gff3', 't.PAC2_0.312.gff3',	't.PAC2_0.316.gff3']

gffFiles = ['t.PAC2_0.189.gff3','t.PAC2_0.296.gff3','t.PAC2_0.311.gff3'] #189 311 296


for gff in gffFiles:
    outFileName = gff[:gff.rfind('.')]+'.sort2'
    inputFile = open(gff, 'r')
    open(outFileName, 'w').close()
    outputFile = open(outFileName, 'w')
    for line in inputFile:
        if 'mRNA' in line and 'longest=1' in line:
            lineInList = line.split()
            parserList = lineInList[-1].split(';')
            lineOutputList = [parserList[1].strip('Name='), lineInList[0], lineInList[3], lineInList[4]]
            outputFile.write('%s    %s   %s   %s\n' % tuple(lineOutputList))

    inputFile.close()
    outputFile.close()
"""

outFileNameN='q.PAC4GC.523.sort'
outFileNameK='q.PAC4GC.524.sort'

inputFile = open('t.PAC2_0.383.gff3', 'r')
open(outFileNameN, 'w').close()
open(outFileNameK, 'w').close()
outputFileN = open(outFileNameN, 'w')
outputFileK = open(outFileNameK, 'w')
for line in inputFile:
    if 'mRNA' in line and 'longest=1' in line:
        lineInList = line.split()
        parserList = lineInList[-1].split(';')
        lineOutputList = [parserList[1].strip('Name='), lineInList[0], lineInList[3], lineInList[4]]
        if 'N' in lineInList[0]:
            outputFileN.write('%s    %s   %s   %s\n' % tuple(lineOutputList))
        if 'K' in lineInList[0]:
            outputFileK.write('%s    %s   %s   %s\n' % tuple(lineOutputList))
        if 'J' in lineInList[0]:
            outputFileN.write('%s    %s   %s   %s\n' % tuple(lineOutputList))
            outputFileK.write('%s    %s   %s   %s\n' % tuple(lineOutputList))

inputFile.close()
outputFileK.close()
outputFileN.close() """