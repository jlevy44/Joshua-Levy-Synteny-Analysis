import subprocess, shutil
"""465-chr2	12171318	13150530	167-Chr5-413519-648466
465-chr3	7744836	8614721	167-Chr5-23212689-23406351
465-chr5	17062130	18365030	167-Chr1-6918190-7208998
465-chr2	12412321	13134727	167-Chr2-15639938-15770619
465-chr1	9950253	10715587	167-Chr2-16754938-17002468
465-chr3	26362082	26902181	167-Chr5-22639843-22792161  """
# convert each entry in above bed file into

# convert to 3 files
# FILE 1 Query
"""Chr05N   84236768	90242904    523_1   100
etc""" # Fakegene from query

# FILE 2 Target
"""Chr05   7419474	10610056    308_1   100
etc""" # fake gene from target

# FILE 3 Mapping the query and target together
"""523_1    308_1   100
etc"""

# files 1-3 are generated from bed file outputs of the synteny analysis, takes large syntenic segments and turns them into
# fake gene bed lines


def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list  files etc
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or info
            configFile.seek(0)
            return line.split()[-1].strip('\n')


open('run.sh','w').close()

# items to be written to shell script
ShellTextList = []

# parse config file
configFile = open('allMAPConfig.txt','r')

# system path to find allmaps installation files
path = parseConfigFindPath('systemPath',configFile)
pythonPath = parseConfigFindPath('pythonPath',configFile) # path of python modules
bedSyntenies = parseConfigFindList('BedFileList', configFile) # path bed files from syntenies
fastaInputFilename = parseConfigFindPath('fastaInputName', configFile) # input genome
fastaOutputName = parseConfigFindPath('fastaOutputName',configFile) # do not include .fasta, will be final output genome
allMAPImageOutputPath = parseConfigFindPath('allMAPImageOutputPath',configFile) # image path of output images
try:
    bedPath = parseConfigFindPath('bedPath',configFile)
except:
    bedPath = ''

# first part of shell bash script
ShellTextList.append("""#!/bin/bash

export PATH=%s:$PATH
export PYTHONPATH=%s:$PYTHONPATH

### creating the input that allmaps needs:
"""%(path,pythonPath))
bedShellList = []

# for each synteny bed file mapping generated from unouts in Synteny analysis
fakeGeneCount = 1
for i in range(len(bedSyntenies)):
    # query and target species protyome Ids
    qspeciesName = bedSyntenies[i][bedSyntenies[i].find('.')+1:bedSyntenies[i].find('-')]
    tspeciesName = bedSyntenies[i][bedSyntenies[i].find('_')+1:bedSyntenies[i].rfind('_')].strip('0.')

    # output files, target and query fakegenes bed file, and the mappings file between the two
    outFileNames = ['%sfakegene.bed%d.txt'%(qspeciesName,i+1),'%sfakegene.bed%d.txt'%(tspeciesName,i+1),
                    '%s_%s.synt.genemap%d.txt'%(qspeciesName,tspeciesName,i+1)]

    for filename in outFileNames:
        open(filename,'w').close()

    # open these output files
    outputFiles=[open(outFileNames[0],'w'),open(outFileNames[1],'w'),
                 open(outFileNames[2],'w')]
    # input bed file, will convert each bed line (a syntenic segment including multiple genes) into a fake gene and mapping between them
    inputFile = open(bedPath+bedSyntenies[i],'r')

    # grab information from the input synteny bed file, and use it to write fake gene bed (2 of them) and mapping bed
    for line in inputFile:
        lineList = line.split()
        lineList2 = lineList[-1].split('-')
        lineList0 = lineList[0].split('-')


        if int(lineList[2])>int(lineList[1]) and int(lineList2[3])>int(lineList2[2]):
            # query output
            outputFiles[0].write('%s\t%s\t%s\t%s_%d\t100\t+\n'
                                 %(lineList0[1].replace('~','S'),lineList[1],lineList[2],lineList0[0],fakeGeneCount))

            # target output
            outputFiles[1].write('%s\t%s\t%s\t%s_%d\t100\t+\n'
                                 % (lineList2[1].replace('~','S'),lineList2[2],lineList2[3],lineList2[0],fakeGeneCount))

            # mapping file write
            outputFiles[2].write('%s_%d\t%s_%d\t100\n'%(lineList2[0],fakeGeneCount,lineList0[0],fakeGeneCount))

            # next gene
            fakeGeneCount += 1

    for file in outputFiles:
        file.close()

    inputFile.close()

    # use fakegene map bed file and individual species bed files and merge them to create another bed file detailing mapping further
    ShellTextList.append("python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --sbed=%s --qbed=%s -o "
                         "fake%s_%s.fileout10000.bed\n"%(outFileNames[2],outFileNames[0],outFileNames[1],qspeciesName,
                                                       tspeciesName))
    # bed shell list adds this detailed mapping
    bedShellList.append('fake%s_%s.fileout10000.bed '%(qspeciesName,tspeciesName))
    #python -m jcvi.assembly.syntenypath bed 473_283.synt.genemap2.txt --switch --scale=10000 --sbed=473fakegene.bed2.txt --qbed=283fakegene.bed2.txt -o fake473_283.fileout10000.bed

# merge multiple maps together from different species
ShellTextList.append('### merge multiple maps:')
ShellTextList.append('\npython -m jcvi.assembly.allmaps mergebed %s-o %s.bed'%(''.join(bedPair for bedPair in
                                                                                        bedShellList ),fastaOutputName))
# run the merged maps together to reconstruct the genome
ShellTextList.append('\n#running allmaps:')
ShellTextList.append('\n\npython -m jcvi.assembly.allmaps path %s.bed '
                     '%s' % (fastaOutputName,fastaInputFilename))

# write to run.sh
shellScriptFile = open('run.sh','w')
shellScriptFile.writelines(ShellTextList)#(''.join(shellText for shellText in ShellTextList))
shellScriptFile.close()

# run allmaps to reconstruct genome
subprocess.call(['sh', 'run.sh'])
#subprocess.call(['mv', '*.pdf','ALLMAPS_Chromosome_Images'])


# copy output pdf images to specific path
if allMAPImageOutputPath.endswith('/'):
    allMAPImageOutputPath = allMAPImageOutputPath[:-1]

# copy output pdf images to specific path
try:
    listFiles=str(subprocess.Popen(['ls','%s'%allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')]],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n')
    for file in listFiles:
        if file.endswith('.pdf'):
            shutil.copy(allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')+1]+file,allMAPImageOutputPath)
except:
    print 'Error: Unable to move output files. Move *.pdf to ALLMAPS_Chromosome_Images'