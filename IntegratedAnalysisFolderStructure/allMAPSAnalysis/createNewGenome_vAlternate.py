import subprocess, shutil

# convert to 3 files
# FILE 1 Query
"""Chr05N   84236768	90242904    523_1   100
etc"""

# FILE 2 Target
"""Chr05   7419474	10610056    308_1   100
etc"""

# FILE 3 Mapping
"""523_1    308_1   100
etc"""

# files 1 and 2 come from sort2 files

# file 3 comes from unout and info is found from gff files

def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of fai files or unOut filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list files/info
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname/info
            configFile.seek(0)
            return line.split()[-1].strip('\n')



open('run.sh','w').close()

# info to be passed to command line
ShellTextList = []

# parse config file
configFile = open('allMAPConfig2.txt','r') #generate in integrated analysis

# system path to find allmaps installation files
path = parseConfigFindPath('systemPath',configFile)
# path of installed modules
pythonPath = parseConfigFindPath('pythonPath',configFile)
unOutPath = parseConfigFindPath('UnOutPath',configFile)
unOutSyntenies = filter(None,str(subprocess.Popen(['ls','%s'%unOutPath],stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
sortPath = parseConfigFindPath('SortPath',configFile)
sortFiles = filter(None,str(subprocess.Popen(['ls', '%s' % sortPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n'))
fastaInputFilename = parseConfigFindPath('fastaInputName', configFile) #only need one fasta to perform analysis
fastaOutputName = parseConfigFindPath('fastaOutputName',configFile) # do not include .fasta at end of this name
allMAPImageOutputPath = parseConfigFindPath('allMAPImageOutputPath',configFile) # output images will be copied here
lociThreshold = int(parseConfigFindPath('loci_threshold',configFile))

# first part of shell bash script
ShellTextList.append("""#!/bin/bash

export PATH=%s:$PATH
export PYTHONPATH=%s:$PYTHONPATH

### creating the input that allmaps needs:
"""%(path,pythonPath))
bedShellList = [] # bed files to be included as output to shell script

def sort2tobed(sort2fname,sortPath):
    """Convert sort2 files to bed.txt files"""
    inputFile = open(sortPath+sort2fname,'r')
    outputName = sort2fname.replace('.sort2','gene.bed.txt').replace('PAC4GC.','').replace('PAC2_0.','').replace('t.','').replace('q.','')
    open(outputName, 'w').close()
    outputFile = open(outputName,'w')
    for line in inputFile:
        if line:
            lineList = line.split()
            # change each line of sort 2 file into line of bed file (slightly different formatting)
            outputFile.write('%s\t%d\t%s\t%s\t100\t+\n'%(lineList[1],int(lineList[2])-1,lineList[3].strip('\n'),lineList[0]))
    outputFile.close()
    inputFile.close()
    return outputName

def sort2todict(sort2fname,sortPath):
    """Change sort2 into a dictionary object with keys being genes of species sort2, allows for ease of access of data"""
    inputFile = open(sortPath+sort2fname,'r')
    dictConverted = {}
    for line in inputFile:
        if line:
            lineList = line.split()
            """dictionary w/ key gene --> (chr,xi,xf)"""
            dictConverted[lineList[0]]=(lineList[1],str(int(lineList[2])-1),lineList[3].strip('\n'))
    inputFile.close()
    return dictConverted

bedFiles = []
"""
for filename in sortFiles:
    if filename:
        bedFiles.append(sort2tobed(filename,sortPath))
        """
dictOfGenes = {}
for file in sortFiles:
    # create dictionary of species, each species entry has their own sort file converted to dictionary for ease access
    dictOfGenes[file.replace('q.','').replace('t.','').split('.')[1]] = sort2todict(file,sortPath)

geneCount = {}
for i in range(len(unOutSyntenies)):
    # find query and target species names from unout files
    qspeciesName = unOutSyntenies[i][unOutSyntenies[i].find('.')+1:unOutSyntenies[i].find('-')]
    tspeciesName = unOutSyntenies[i][unOutSyntenies[i].find('_')+1:unOutSyntenies[i].rfind('_')].strip('0.')
    # initialize genecounts for each species if doesn't already exist. If exists, continue gene count
    if qspeciesName not in geneCount.keys():
        geneCount[qspeciesName]=0
    if tspeciesName not in geneCount.keys():
        geneCount[tspeciesName]=0
    # replaces sort2bed file, gnerate bed file for each sort file, multiple bed files for query species
    for file in sortFiles:
        if qspeciesName in file:
            qbedFname = file.replace('.sort2','gene.bed%d.txt'%(i+1)).replace('PAC4GC.','').replace('PAC2_0.','').replace('t.','').replace('q.','')
        if tspeciesName in file:
            tbedFname = file.replace('.sort2','gene.bed.txt').replace('PAC4GC.','').replace('PAC2_0.','').replace('t.','').replace('q.','')

    # output file names
    outFileNames = ['%s'%(qbedFname),'%s'%(tbedFname),
                    '%s_%s.synt.genemap%d.txt'%(qspeciesName,tspeciesName,i+1)]

    # create files for writing
    for filename in outFileNames:
        open(filename,'w').close()

    #establish mapping between output files 1 and 2 by generating file 3 from unout
    # unout maps the relationships between species 1 and 2, unout info is in 3, the genemap
    outputFileGenemap=open(outFileNames[2],'w')
    inputFile = open(unOutPath+unOutSyntenies[i],'r')

    read=0
    # output query and target file objects
    outq= open(outFileNames[0],'w')
    outt= open(outFileNames[1],'w')
    for line in inputFile:
        if line:
            if '[' not in line and read == 1: #	Chr05N,3014546,Pavir.5NG471200.1,Chr05,1013018,Pahal.E01839.1	3:1	5	115	782	0.0	0.871794871794872
                # (target species gene name, query species gene name)
                genes = (line.split()[0].split(',')[-1],line.split()[0].split(',')[2])
                # add to genecounts to establish mapping, will rename gene names to reflect gene count but still grab sequence information
                geneCount[tspeciesName] += 1
                geneCount[qspeciesName] += 1
                # the genes found in unout have been renamed into "fakegenes", will still carry same sequence info
                # (target species fakegene name, query species fakegene name)
                fakegenes = ('%s_%d' % (tspeciesName, geneCount[tspeciesName]), '%s_%d' % (qspeciesName, geneCount[qspeciesName]))
                # write that fake gene has information of actual gene
                outq.write('%s\t%s\t%s\t%s\t100\t+\n'%(dictOfGenes[qspeciesName][genes[1]]+(fakegenes[1],)))
                outt.write('%s\t%s\t%s\t%s\t100\t+\n'%(dictOfGenes[tspeciesName][genes[0]]+(fakegenes[0],)))
                #print qspeciesName,dictOfGenes[qspeciesName].has_key(genes[1]),tspeciesName,dictOfGenes[tspeciesName].has_key(genes[0])
                # establish relationship between two syntenic genes through genemap file
                outputFileGenemap.write('%s\t%s\t100\n'%fakegenes)
            if '[' in line:
                read =0
                # only take in genes that pass loci threshold
                if int(line.split()[-6])>=lociThreshold:
                    read =1


    outputFileGenemap.close()
    outq.close()
    outt.close()

    inputFile.close()

    # use gene map bed file and individual species bed files and merge them to create another bed file detailing mapping further
    ShellTextList.append("python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --sbed=%s --qbed=%s -o "
                         "%s_%s.fileout10000.bed\n"%(outFileNames[2],outFileNames[0],outFileNames[1],qspeciesName,
                                                       tspeciesName))
    # bed shell list adds this detailed mapping
    bedShellList.append('%s_%s.fileout10000.bed '%(qspeciesName,tspeciesName))
    #python -m jcvi.assembly.syntenypath unOut 473_283.synt.genemap2.txt --switch --scale=10000 --sunOut=473fakegene.unOut2.txt --qunOut=283fakegene.unOut2.txt -o fake473_283.fileout10000.unOut

# merge multiple maps together from different species
ShellTextList.append('### merge multiple maps:')
ShellTextList.append('\npython -m jcvi.assembly.allmaps mergebed %s-o %s.bed'%(''.join(bedPair for bedPair in
                                                                                        bedShellList ),fastaOutputName))
# run the merged maps together to reconstruct the genome
ShellTextList.append('\n#running allmaps:')
ShellTextList.append('\n\npython -m jcvi.assembly.allmaps path %s.bed '
                     '%s' % (fastaOutputName,fastaInputFilename))

# write shell text list to run.sh
shellScriptFile = open('run.sh','w')
shellScriptFile.writelines(ShellTextList)#(''.join(shellText for shellText in ShellTextList))
shellScriptFile.close()

# run allmaps to reconstruct genome
subprocess.call(['sh', 'run.sh'])
#subprocess.call(['mv', '*.pdf','ALLMAPS_Chromosome_Images'])


# copy output pdf images to specific path
if allMAPImageOutputPath.endswith('/'):
    allMAPImageOutputPath = allMAPImageOutputPath[:-1]

# copy images
try:
    listFiles=str(subprocess.Popen(['ls','%s'%allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')]],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n')
    for file in listFiles:
        if file.endswith('.pdf'):
            shutil.copy(allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')+1]+file,allMAPImageOutputPath)
except:
    print 'Error: Unable to move output files. Move *.pdf to ALLMAPS_Chromosome_Images'

