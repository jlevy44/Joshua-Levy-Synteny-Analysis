import subprocess, shutil
"""465-chr2	12171318	13150530	167-Chr5-413519-648466
465-chr3	7744836	8614721	167-Chr5-23212689-23406351
465-chr5	17062130	18365030	167-Chr1-6918190-7208998
465-chr2	12412321	13134727	167-Chr2-15639938-15770619
465-chr1	9950253	10715587	167-Chr2-16754938-17002468
465-chr3	26362082	26902181	167-Chr5-22639843-22792161  """

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


def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of fai files or bed filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of fai or bed files
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname
            configFile.seek(0)
            return line.split()[-1].strip('\n')


open('run.sh','w').close()

ShellTextList = []

# parse config file
configFile = open('allMAPConfig.txt','r')

path = parseConfigFindPath('systemPath',configFile)
pythonPath = parseConfigFindPath('pythonPath',configFile)
bedSyntenies = parseConfigFindList('BedFileList', configFile)
fastaInputFilename = parseConfigFindPath('fastaInputName', configFile)
fastaOutputName = parseConfigFindPath('fastaOutputName',configFile)
allMAPImageOutputPath = parseConfigFindPath('allMAPImageOutputPath',configFile)
try:
    bedPath = parseConfigFindPath('bedPath',configFile)
except:
    bedPath = ''

ShellTextList.append("""#!/bin/bash

export PATH=%s:$PATH
export PYTHONPATH=%s:$PYTHONPATH

### creating the input that allmaps needs:
"""%(path,pythonPath))
bedShellList = []

fakeGeneCount = 1
for i in range(len(bedSyntenies)):
    qspeciesName = bedSyntenies[i][bedSyntenies[i].find('.')+1:bedSyntenies[i].find('-')]
    tspeciesName = bedSyntenies[i][bedSyntenies[i].find('_')+1:bedSyntenies[i].rfind('_')].strip('0.')

    outFileNames = ['%sfakegene.bed%d.txt'%(qspeciesName,i+1),'%sfakegene.bed%d.txt'%(tspeciesName,i+1),
                    '%s_%s.synt.genemap%d.txt'%(qspeciesName,tspeciesName,i+1)]

    for filename in outFileNames:
        open(filename,'w').close()


    outputFiles=[open(outFileNames[0],'w'),open(outFileNames[1],'w'),
                 open(outFileNames[2],'w')]
    inputFile = open(bedPath+bedSyntenies[i],'r')



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

    ShellTextList.append("python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --sbed=%s --qbed=%s -o "
                         "fake%s_%s.fileout10000.bed\n"%(outFileNames[2],outFileNames[0],outFileNames[1],qspeciesName,
                                                       tspeciesName))
    bedShellList.append('fake%s_%s.fileout10000.bed '%(qspeciesName,tspeciesName))
    #python -m jcvi.assembly.syntenypath bed 473_283.synt.genemap2.txt --switch --scale=10000 --sbed=473fakegene.bed2.txt --qbed=283fakegene.bed2.txt -o fake473_283.fileout10000.bed

ShellTextList.append('### merge multiple maps:')
ShellTextList.append('\npython -m jcvi.assembly.allmaps mergebed %s-o %s.bed'%(''.join(bedPair for bedPair in
                                                                                        bedShellList ),fastaOutputName))
ShellTextList.append('\n#running allmaps:')
ShellTextList.append('\n\npython -m jcvi.assembly.allmaps path %s.bed '
                     '%s' % (fastaOutputName,fastaInputFilename))

# write to run.sh
shellScriptFile = open('run.sh','w')
shellScriptFile.writelines(ShellTextList)#(''.join(shellText for shellText in ShellTextList))
shellScriptFile.close()

subprocess.call(['sh', 'run.sh'])
#subprocess.call(['mv', '*.pdf','ALLMAPS_Chromosome_Images'])

if allMAPImageOutputPath.endswith('/'):
    allMAPImageOutputPath = allMAPImageOutputPath[:-1]

try:
    listFiles=str(subprocess.Popen(['ls','%s'%allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')]],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n')
    for file in listFiles:
        if file.endswith('.pdf'):
            shutil.copy(allMAPImageOutputPath[:allMAPImageOutputPath.rfind('/')+1]+file,allMAPImageOutputPath)
except:
    print 'Error: Unable to move output files. Move *.pdf to ALLMAPS_Chromosome_Images'