from bed2linkfile import *
from generateConfigs import *
from fai2karyotype import *
import subprocess
import sys


# DATA
# following functions parse the configuration file
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
    configFile.seek(0)


# time to parse config file
configFile = open('circosFiguresPipelineConfig.txt','r')
for line in configFile:
    if 'pathPythonModules' in line: # add python path in case not in root directory, specifies where to find modules/lib
        sys.path.append(line.split()[-1].strip('\n'))
        configFile.seek(0)
        break

faiFiles = parseConfigFindList('faiFiles',configFile) # find fai files as list, will turn into karyotypes
listBedFiles = parseConfigFindList('BedFiles',configFile) # bed file list, will turn into link files
findPathsList = ['faiFilePath','karyotypesFilesPath','BedInputPath','circosConfigFilesPath','LinkPath',
                 'circosOutPath','BPsThreshold'] # list of query strings into config path finder
for i in range(len(findPathsList)): # find the paths of above queries
    findPathsList[i] = parseConfigFindPath(findPathsList[i],configFile)
# associate pathnames with individual variables, BPs threshold is number of basepairs to accept as cutoff for adding
# chromosome to karyotype file
(faiFilePath, karyotypesFilesPath, BedInputPath, circosConfigFilesPath, LinkPath,circosOutPath, BPsThreshold)=\
    tuple(findPathsList)


# convert bed files into link files
linkFiles = bed2link(listBedFiles,BedInputPath,LinkPath)

# find name of query species eg. 283 # FIXME this kind of does not work for comparisons for 2 or more query species
querySpecies = linkFiles[0][linkFiles[0].find('.')+1:linkFiles[0].find('-')]


# find fai filename associated with query species
for faiFile in faiFiles:
    if querySpecies in faiFile:
        mainSpeciesFai = faiFile

# convert fai files to karyotype files
listKaryotypes = fai2karyotype(faiFilePath,karyotypesFilesPath,faiFiles,mainSpeciesFai,BPsThreshold)

# for each link file, generate a circos picture describing pairwise comparison synteny between species
for linkFile in linkFiles:
    KaryotypeFiles = []
    speciesList = [0,0]
    # name of query and target species
    speciesList[0] = linkFile[linkFile.find('.')+1:linkFile.find('-')]
    speciesList[1] = linkFile[linkFile.find('-')+1:linkFile.rfind('_')].replace('PAC2_0.','').replace('PAC4GC.','')
    # add 2 karyotypes file names to list describing target and query species for pairwise comparison
    for species in speciesList:
        for karyotype in listKaryotypes:
            if species in karyotype:
                KaryotypeFiles.append(karyotype)
    # P virgatum (can modify later) has 2 subgenomes, take this into account in final analysis
    if '523' in speciesList and '524' in speciesList:
        subgenomesCheck = 1
    else:
        subgenomesCheck = 0
    # generate configuration files for circos plot
    generateConfigs(circosConfigFilesPath,karyotypesFilesPath,LinkPath,KaryotypeFiles,linkFile,subgenomesCheck)
    # create circos plot from the command line and export to specified path
    print circosConfigFilesPath, speciesList, circosOutPath
    subprocess.call(['circos','-conf',circosConfigFilesPath+'circos.conf','-outputfile',
                     '%s-%s'%(speciesList[0],speciesList[1]),'-outputdir',circosOutPath])
    #subprocess.call('circos -conf %scircos.conf -outputfile %s-%s -outputdir %s'%(circosConfigFilesPath,speciesList[0],
     #                                                                             speciesList[1],circosOutPath),shell=True)
    #print '/Applications/circos-0.69/bin/circos -conf %s'%(configPath+'circos.conf')
    #/Applications/circos-0.69/bin/