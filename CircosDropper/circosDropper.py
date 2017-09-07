from fai2karyotypeMod import fai2karyotype
from bed2linkfileMod import bed2link
from SyntenyFinalMod import unout2bed
from generateConfigsMod import generateConfigs
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
import os, subprocess
def gff2sort2(gff, pathgff='./', pathsort='./'):
    """Takes a gffFiles and converts them to sort2 files to use in the final synteny analysis.
    Please let Joshua Levy know if there are any errors or have problems!"""
    outFileName = pathsort + gff[:gff.rfind('.')] + '.sort2'
    inputFile = open(pathgff + gff, 'r')
    open(outFileName, 'w').close()
    outputFile = open(outFileName, 'w')
    for line in inputFile:
        # grab gene info from each line if it's longest and mRNA strand and output to sort2 file
        if 'mRNA' in line and 'longest=1' in line:
            lineInList = line.split()
            parserList = lineInList[-1].split(';')
            lineOutputList = [parserList[1].replace('Name=',''), lineInList[0].replace('-', 'S'), lineInList[3],
                              lineInList[4]]
            outputFile.write('%s    %s   %s   %s\n' % tuple(lineOutputList))

    inputFile.close()
    outputFile.close()

findFiles = lambda x: [file for file in os.listdir('.') if file.endswith(x)]

gffFiles = findFiles('.gff3')

sortFiles = [file.replace('.gff3','.sort2') for file in gffFiles]
for sortFile in sortFiles:
    if sortFile not in os.listdir('.'):
        gff2sort2(sortFile.replace('.sort2','.gff3'))


karyotypeFiles =  [fai2karyotype(faifile) for faifile in findFiles('.fai')]

unoutFiles = findFiles('.unout')

for unoutfile in unoutFiles:
    bedfile = unoutfile.replace('.unout','.bed')
    querySpecies = unoutfile[unoutfile.find('.')+1:unoutfile.find('-')]
    targetSpecies = unoutfile[unoutfile.find('-') + 1:unoutfile.rfind('_')].replace('PAC2_0.', '').replace('PAC4GC.', '')
    #print targetSpecies
    for sortFile in sortFiles:
        if querySpecies in sortFile:
            querySortFile = sortFile
        if targetSpecies in sortFile:
            targetSortFile = sortFile
            #print targetSortFile
    if bedfile not in os.listdir('.'):
        try:
            unout2bed((unoutfile,querySortFile,targetSortFile))
        except:
            print (unoutfile,querySortFile,targetSortFile)
            quit()
    linkFile = unoutfile.replace('.unout','.link.txt')
    if linkFile not in os.listdir('.'):
        bed2link(bedfile)
    karyotypeQT = tuple([karyFile for karyFile in karyotypeFiles if querySpecies in karyFile.split('.')[1] or targetSpecies in karyFile.split('.')[1]])
    generateConfigs(karyotypeQT,linkFile)
    if '%s-%s.png'%(querySpecies,targetSpecies) not in os.listdir('.'):
        #try:
        subprocess.call('circos -conf circos.conf -outputfile '+'%s-%s' %(querySpecies,targetSpecies)+ ' -outputdir .',shell=True)
        #except:
        #    print 'circos', '-conf', 'circos.conf', '-outputfile', '%s-%s' %(querySpecies,targetSpecies), '-outputdir', '.'
        #    exit()
    if '%s-%s.pdf' % (querySpecies, targetSpecies) not in os.listdir('.'):
        subprocess.call('convert %s-%s.png %s-%s.pdf'%(querySpecies,targetSpecies,querySpecies,targetSpecies),shell=True)
        #try:
        #    drawing = svg2rlg('%s-%s.svg'%(querySpecies,targetSpecies))
        #    renderPDF.drawToFile(drawing,'%s-%s.pdf' % (querySpecies, targetSpecies))
        #except:
        #    print '%s-%s.pdf' % (querySpecies, targetSpecies)