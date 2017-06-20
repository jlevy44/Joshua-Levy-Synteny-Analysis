from generateIntegratedConfigurationFiles import generateIntAnalysisConfig
import subprocess, os, shutil, sys
from multiprocessing import Pool
"""
try:
    import syntenyFinal
except:
    print'cannot find syntenyFinal analysis'
    exit()
try:
    import circosFiguresPipeline
except:
    print'cannot find circos analysis'
    exit()
    """
# gff2sort before completing entire analysis
# name bed files folder BdFiles


# PARSING FUNCTIONS

def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of relevant filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of files or list information
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    configFile.seek(0)
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path or value of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or specific value trying to find
            configFile.seek(0)
            return line.split()[-1].strip('\n')
    configFile.seek(0)


# open master configuration file
open('masterConfig.txt','r').close()
masterConfigFile = open('masterConfig.txt','r')


# grab the following information from the configuration file
weightsList = parseConfigFindList('Weights info',masterConfigFile)
findInfoList = ['performSynteny','performCircos', 'performALLMAPS', 'querySpecies', 'NameAnalysis','writeFastaOut', 'Loci_Threshold',
                'pathPython','pathSystem', 'pathALLMAPS', 'BdPath', 'pathUnOut', 'pathGFF', 'pathSort', 'BPsMergeDist', 'softMasked', 'genomePath',
                'karyotypesFilesPath','circosConfigFilesPath', 'LinkPath', 'circosOutPath', 'BPsThreshold',
                'multipleSeqAlignFastasPath','fastaOutputName', 'allMAPImageOutputPath', 'online','projectName',
                'nerscUsername','nohup','cactusRun','cactusFolder'] # find the following information
 # list of query strings into config path finder
for i in range(len(findInfoList)): # find the paths/info of above queries
    findInfoList[i] = parseConfigFindPath(findInfoList[i], masterConfigFile)

# assign values
(performSynteny, performCircos, performALLMAPS,querySpecies, NameAnalysis, writeFastaOut, Loci_Threshold, pathPython,
 pathSystem,pathALLMAPS, BdPath, pathUnOut, pathGFF, pathSort, BPsMergeDist , softmask, genomePath, karyotypesFilesPath,
 circosConfigFilesPath, LinkPath, circosOutPath, BPsThreshold, multipleSeqAlignFastasPath,
 fastaOutputName, allMAPImageOutputPath, online, projectName, nerscUsername, nohup, cactusRun,cactusFolder) = tuple(findInfoList)



# for debugging, see if all of your data has passed through
print tuple(findInfoList)

# generate weights file for allmaps...
open('%sweights.txt'%pathALLMAPS,'w').close()
weightsFile = open('%sweights.txt'%pathALLMAPS,'w')
for weight in weightsList:
    weightsFile.write('fake%s_'%querySpecies+weight+'\n')
weightsFile.close()
#fake473_283 1

# generate config files
# first need to generate 3 syntenic files list and list of genome files for synteny analysis
# second generate faiFiles list and bed files list circos
# third generate bedfile list and fastainputname

# get list of genome files, fai files, and fasta input filename for fragmented genome
listGenomePathFiles=str(subprocess.Popen(['ls','%s'%genomePath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        .stdout.read()).split('\n')
listGenomeFiles = []
listFaiFiles = []

# find lists of .fa and .fai files and format for config files
for file in listGenomePathFiles:
    if file.endswith('.fa') or file.endswith('.fasta'):
        listGenomeFiles.append(file+'\n')
        listFaiFiles.append(file+'.fai\n')

genomeFilesText = ''.join(file for file in listGenomeFiles)
faiFilesText = ''.join(file for file in listFaiFiles)

# if query species, use .fa file for genome reconstruction
for filename in listGenomeFiles:
    if querySpecies in filename:
        fastaInputName = filename.strip('\n') # fasta filename of query species

# list of unout files
listPathUnout = str(subprocess.Popen(['ls', '%s' % pathUnOut], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')


def gff2sort2(gff, pathgff, pathsort):
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

# prints True if there are no sortfiles, generally want that for each analysis... for now, is possible to modify code to
# be more versatile
print not str(subprocess.Popen(['ls', '%s' % pathSort], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read())

# if no sort files present, generate them from gffs. PLEASE DELETE SORT FILES AFTER EVERY ANALYSIS
if not str(subprocess.Popen(['ls', '%s' % pathSort], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()):

    listGFFfiles = str(subprocess.Popen(['ls', '%s' % pathGFF], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                       .stdout.read()).split('\n')
    # turn gff into sort files
    for file in listGFFfiles:
        if file.endswith('.gff') or file.endswith('.gff3'):
            gff2sort2(file,pathGFF,pathSort)

# find sort files
listPathSort = str(subprocess.Popen(['ls', '%s' % pathSort], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                          .stdout.read()).split('\n')

unoutList = []
sortFileList = []
bedList = []

# generate list of bed files for allmaps and circos, and unout and sort files for synteny analysis
# unouts can also be used alternate allmaps
for file in listPathUnout:
    if file.endswith('.unout'):
        bedList.append(file[:file.rfind('.')]+'.bed')
        unoutList.append(file)
for file in listPathSort:
    if file.endswith('.sort2'):
        sortFileList.append(file)

# make sure all sort files are included
print sortFileList

# Bedfile text for config files
bedFilesText = ''.join(file+'\n' for file in bedList)

# generate tuples for three syntenic files
listSyntenicFilesTuples = []
for file in unoutList:
    # find target species then corresponding target sort files, add to list of syntenic files
    targetSpecies = file[file.find('-')+1:file.rfind('_')].replace('PAC2_0.','').replace('PAC4GC.','')
    print targetSpecies
    for sortFile in sortFileList:
        if querySpecies in sortFile:
            querySortFile = sortFile
        if targetSpecies in sortFile:
            targetSortFile = sortFile
            print targetSortFile
    listSyntenicFilesTuples.append((file,querySortFile,targetSortFile))

# text for synteny analysis config
syntenicFilesText = ''.join('%s   %s   %s\n'%synFilesTuple for synFilesTuple in listSyntenicFilesTuples)

# if on NERSC load proper modules, try to load these modules beforehand...
print 'online' + online
if int(online):
    try:
        subprocess.call('module load bedtools/2.25.0',shell=True)
        subprocess.call('module load circos',shell=True)
    except:
        print 'Unable to load online modules...'

try:
    int(BPsMergeDist)
except:
    BPsMergeDist = '100000'

# write syntenic Config text
generateIntAnalysisConfig('syntenyAnalysis',(NameAnalysis,writeFastaOut,Loci_Threshold,pathPython,pathUnOut,pathSort, BPsMergeDist,softmask,
                                             NameAnalysis,multipleSeqAlignFastasPath,syntenicFilesText,genomePath,
                                             genomeFilesText))

# if choosing to perform multiple synteny
if int(performSynteny):
    # run synteny analysis
    #execfile(os.path.join(os.path.dirname(sys.argv[0]), 'syntenyFinal.py'))
    execfile('SyntenyFinal.py')
    # move BedFiles to Bd folder
    for bedFile in bedList:
        try:
            shutil.copy(BdPath[:BdPath.rfind('BdFiles')]
                    + bedFile, BdPath)
        except:
            print 'Unable to copy bed file %s to BdFiles...'%bedFile
    # move genome to allmaps folder
    try:
        shutil.copy(genomePath+fastaInputName,pathALLMAPS+fastaInputName)
    except:
        print 'Unable to copy genome to allmaps directory...'


# generate config files for relevant analyses
if int(performCircos):
    print 'config circos'
    generateIntAnalysisConfig('circosAnalysis',(faiFilesText,bedFilesText,genomePath,karyotypesFilesPath,BdPath,
                                                circosConfigFilesPath,LinkPath,circosOutPath,pathPython,BPsThreshold))
if int(performALLMAPS):
    print 'config allmaps'
    generateIntAnalysisConfig('reconstructGenome',(pathALLMAPS,pathPython,pathSystem,BdPath,allMAPImageOutputPath,
                                                   fastaInputName,fastaOutputName,bedFilesText))

# see if will be performing alternate genome reconstruction using syntenic genes, MORE ACCURATE
try:
    masterConfigFile = open('masterConfig.txt','r')
    masterConfigFile.seek(0)
    performAltAllMaps = parseConfigFindPath('performAltAllMaps', masterConfigFile)
    masterConfigFile.close()
    if int(performAltAllMaps):
        print 'config alt allmaps'
        generateIntAnalysisConfig('reconstructGenome2', (pathALLMAPS, pathPython, pathSystem, BdPath, allMAPImageOutputPath,
                                                fastaInputName, fastaOutputName, pathUnOut,pathSort,Loci_Threshold))
        open('%sweights.txt' % pathALLMAPS, 'w').close()
        weightsFile = open('%sweights.txt' % pathALLMAPS, 'w')
        for weight in weightsList:
            weightsFile.write('%s_' % querySpecies + weight + '\n')
        weightsFile.close()
except:
    print 'Unable to set up alternate allmaps'
    performAltAllMaps = 0

# perform circos analysis but not reconstruction
if int(performCircos) and not int(performALLMAPS):
    print 'circos'
    # try to run circos online
    if int(online):
        open('runCircos.sh', 'w').close()
        # writing shell script to run circos
        circos = open('runCircos.sh', 'w')
        circos.write('#!/bin/bash\npython circosFiguresPipeline.py')
        circos.close()
        try:
            # try to run circos
            subprocess.call('nohup sh runCircos.sh', shell=True)
        except:
            print 'Unable to run circos via command line..'
    else: # try offline analysis
        try:
            execfile('circosFiguresPipeline.py')
        except:
            print 'Unable to run circos analysis.'
    #except:
    #    print 'Unable to run circos analysis.'
    #    exit()

# perform allmaps but not circos
elif int(performALLMAPS) and not int(performCircos):
    print 'allmaps'
    #os.chdir(pathALLMAPS)
    # if online
    if int(online):
        # try to run allmaps
        allmap = open('runAllmaps.sh', 'w')
        allmap.write( '#!/bin/bash\ncd %s\npython createNewGenome.py'%pathALLMAPS)
        allmap.close()
        try:
            if int(nohup): # if using no hangup on a server
                subprocess.call('nohup sh runAllmaps.sh', shell=True)
            else: # else submit the job to run on the cluster
                subprocess.call('qsub -P %s -N allmapsAnalyze -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
                 ' -l exclusive.c %s%s\n'%(projectName,BdPath[:BdPath.rfind('BdFiles')],nerscUsername,'runAllmaps.sh'),
                            shell=True)
        except:
            print 'Unable to run allmaps via qsub. Set online to 0 and try again. Trying to run via command line..'
            try:
                subprocess.call('nohup sh runAllmaps.sh', shell=True)
            except:
                print 'Unable to run allmaps via command line..'
    # else run in offline mode
    else:
        try:
            os.chdir(pathALLMAPS)
            execfile('createNewGenome.py')
        except:
            print 'Unable to run allmaps analysis.'
            exit()
# perform both allmaps and circos
elif int(performALLMAPS) and int(performCircos):
    print 'both'

    if int(online):
        open('runCircos.sh','w').close()
        circos = open('runCircos.sh','w')
        circos.write('#!/bin/bash\npython circosFiguresPipeline.py')
        circos.close()
        # multiple command line scripts being written to run on the cluster, talk to NERSC about syntax
        lines = ['qsub -P %s -N circosAnalyze -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
                 ' -l exclusive.c %s%s\n'%(projectName,BdPath[:BdPath.rfind('BdFiles')],nerscUsername,'runCircos.sh'),
                 'cd %s\n'%pathALLMAPS,
                 'qsub -P %s -N allmapsAnalyze -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
                 ' -l exclusive.c %s%s\n'%(projectName,BdPath[:BdPath.rfind('BdFiles')],nerscUsername,'runAllmaps.sh')]
        open('runAllmaps.sh', 'w').close()
        allmap = open('runAllmaps.sh', 'w')
        allmap.write('#!/bin/bash\n'+lines[1]+'python createNewGenome.py')
        allmap.close()

        try: # run circos on nohup
            subprocess.call('nohup sh runCircos.sh', shell=True)
        except:
            print 'Unable to run circos analysis.'
        try:
            #subprocess.call(lines[0],shell=True)
            if int(nohup): # run allmaps on nohup
                subprocess.call('nohup sh runAllmaps.sh', shell=True)
            else: # try submission to cluster if does not work
                subprocess.call(lines[2],shell=True)
        except:
            print 'Unable to run allmaps via qsub. Set online to 0 and try again. Trying to run via command line..'
            try:
                subprocess.call('nohup sh runAllmaps.sh',shell=True)
            except:
                print 'Unable to run allmaps via command line..'
    else: # offline mode, linux preferred
        """
        def runALLMAPS():
            os.chdir(pathALLMAPS)
            try:
                execfile('createNewGenome.py')
            except:
                print 'Unable to run allmaps analysis.'
        pool = Pool(processes=2)
        try:
            pool.apply_async(execfile,'circosFiguresPipeline.py')
        except:
            print 'Unable to run circos analysis.'
        pool.apply_async(runALLMAPS)
        pool.close()
        pool.join()
        """
        try:
            execfile('circosFiguresPipeline.py')
        except:
            print 'Unable to run circos analysis.'
        os.chdir(pathALLMAPS)
        try:
            print 'Warning: If online, monitor computer for no timeout...'
            execfile('createNewGenome.py')
        except:
            print 'Unable to run allmaps analysis.'
            exit()

if int(performAltAllMaps): # if perform alternative ALLMAPS, RECOMMENDED over regular allmaps, similar to allmaps but
    # not circos code
    if int(online):
        allmap = open('runAllmaps2.sh', 'w')
        allmap.write( '#!/bin/bash\n\ncd %s\npython createNewGenome_vAlternate.py'%pathALLMAPS)
        allmap.close()
        try:
            if int(nohup):
                subprocess.call('nohup sh runAllmaps2.sh', shell=True)
            else:
                subprocess.call('qsub -P %s -N allmapsAnalyze -cwd -b yes -now no -j yes -m abes -M %s@lbl.gov -w e'
                 ' -l exclusive.c %s%s\n'%(projectName,BdPath[:BdPath.rfind('BdFiles')],nerscUsername,'runAllmaps2.sh'),
                            shell=True)
        except:
            print 'Unable to run allmaps via qsub. Set online to 0 and try again. Trying to run via command line..'
            try:
                subprocess.call('nohup sh runAllmaps2.sh', shell=True)
            except:
                print 'Unable to run alt allmaps from command line'
    else:
        try:
            os.chdir(pathALLMAPS)
            execfile('createNewGenome_vAlternate.py')
        except:
            print 'Unable to run allmaps analysis.'
            exit()

try:
    if int(cactusRun):
        for file in os.listdir(multipleSeqAlignFastasPath):
            if file.endswith('fasta'):
                shutil.copy(file,cactusFolder+'FastaFiles')
        print 'Remember to have grassnewick file in place!!!'
        os.chdir(cactusFolder)
        subprocess.call('sh runcac.sh',shell = True)
except:
    print 'Error running cactus...'