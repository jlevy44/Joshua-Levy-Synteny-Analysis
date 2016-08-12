from pybedtools import *

bedList = ['PAC4GC.524-PAC2_0.308_5.bed','PAC4GC.524-PAC2_0.312_5.bed','PAC4GC.524-PAC2_0.313_5.bed',
           'PAC4GC.524-PAC2_0.314_5.bed','PAC4GC.524-PAC2_0.316_5.bed','PAC4GC.524-PAC2_0.323_5.bed',
           'PAC4GC.524-PAC4GC.523_5.bed']

finalBedIntersect = BedTool(bedList[0])
for bed in bedList[1:]:
    finalBedIntersect = finalBedIntersect.intersect(BedTool(bed)).sort().merge()

finalBedMerge = BedTool(bedList[0])
for bed in bedList[1:]:
    finalBedMerge = finalBedMerge.cat(BedTool(bed), postmerge = False)

finalBedMerge = finalBedMerge.sort().merge()
print finalBedIntersect

print finalBedMerge

# option to create histogram for pairwise comparison...
#finalBedMerge = BedTool('PAC2_0.283-PAC2_0.323_5.bed').sort().merge()

inputFile = open('GENE_LENGTH_Pvirgatum_383_v3.0.softmasked.filt.K.bed','r')
#open('all.shell.core.TE.low_depth.SNPs.recombRate.TEmatrixINS.TEmatrixABS.synteny.filt.bedgraph.ShellRatio.'
                 #'tr.txt','r')
open('histogramIntervalsBed.bed','w').close()
histogramIntFile = open('histogramIntervalsBed.bed','w')
#inputFileLines = inputFile.readlines()
#for line in inputFileLines[1:]:
#    histogramIntFile.write('523-%s\t%s\t%s\n'%(tuple(line.split()[0:3])))

for line in inputFile:
    histogramIntFile.write('524_'+line)

histogramIntFile.close()
inputFile.close()

histogramBed = BedTool('histogramIntervalsBed.bed').sort()

print histogramBed

#for line in finalBedMerge:
#    print line
#finalHistFile =

outIntersectHist = open('524.syntBP.conservedAmongAllSpecies.histogram.txt','w')
outMergeHist = open('524.syntBP.mergedOverlappingIntervalsAllSpecies.histogram.txt','w')
#mergedOverlappingIntervalsAllSpecies
totalFiles = 2
pairwise = 1


try:
    beds = [histogramBed.intersect(finalBedMerge,wao = True).merge(c=7,o='distinct',d=-3),
            histogramBed.intersect(finalBedIntersect, wao=True).merge(c=7, o='distinct', d=-3)]
except:
    totalFiles = 1
    beds = [histogramBed.intersect(finalBedMerge, wao=True).merge(c=7, o='distinct', d=-3)]
histFiles = [outMergeHist, outIntersectHist]

for i in range(totalFiles):
    for line in beds[i]:
        lineList = str(line).split()
        synBPsList = lineList[-1].split(',')
        synBPsList[-1].strip('/n')
        for j in range(len(synBPsList)):
            synBPsList[j] = int(synBPsList[j])
        histFiles[i].write('%s %s %s %f\n' % (lineList[0][lineList[0].find('-')+1:],lineList[1],lineList[2],
                                              sum(synBPsList)/float(int(lineList[2])-int(lineList[1]))))

for file in histFiles:
    file.close()

#print histogramBed.intersect(finalBedIntersect,wao = True).merge(c=7,o='distinct',d=-3)#.sort().merge()
print histogramBed.intersect(finalBedMerge,wao=True).merge(c=7,o = 'distinct',d=-3)#,wao = True).sort().merge()

"""
for intervalLine in histogramBed:
    intervalList = str(intervalLine).split()
    synBPs = [[],[]]
    totalBPs = int(intervalList[2])-int(intervalList[1])
    for i in range(2):
        for line in beds[i]:
            if intervalLine[0] in line and ( int(intervalList[1])<=int(str(line).split()[1])<=int(intervalList[2])  or  \
                                    int(intervalList[1])<=int(str(line).split()[2])<=int(intervalList[2]) ):

                synBPs[i] += [len(set(range(int(str(line).split()[1]),int(str(line).split()[2])))
                                  .intersection(range(int(intervalList[1]),int(intervalList[2]))))]

        synBPsRatio = sum(synBPs[i])/float(totalBPs)

        histFiles[i].write('%s %s %s %d\n'%(intervalList[0][intervalList[0].find('-')+1:],intervalList[1],intervalList[2],
                                          synBPsRatio))
"""