from bed2linkfile import *
from generateConfigs import *
from fai2karyotype import *
import subprocess
import sys


# DATA
def parseConfigFindList(stringFind,configFile):
    read = 0
    listOfItems = []
    for line in configFile:
        if read == 1:
            if 'Stop' in line:
                configFile.seek(0)
                break
            listOfItems.append(line.strip('\n'))
        if stringFind in line:
            read = 1
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    for line in configFile:
        if stringFind in line:
            return line.split()[-1].strip('\n')


# time to parse config file
configFile = open('circosFiguresPipelineConfig.txt','r')
for line in configFile:
    if 'pathPythonModules' in line:
        sys.path.append(line.split()[-1].strip('\n'))
        configFile.seek(0)
        break
faiFiles = parseConfigFindList('faiFiles',configFile)
listBedFiles = parseConfigFindList('BedFiles',configFile)
findPathsList = ['faiFilePath','karyotypesFilesPath','BedInputPath','circosConfigFilesPath','LinkPath',
                 'circosOutPath']
for i in range(len(findPathsList)):
    findPathsList[i] = parseConfigFindPath(findPathsList[i],configFile)
(faiFilePath, karyotypesFilesPath, BedInputPath, circosConfigFilesPath, LinkPath,circosOutPath)=\
    tuple(findPathsList)



linkFiles = bed2link(listBedFiles,BedInputPath,LinkPath)

querySpecies = linkFiles[0][linkFiles[0].find('.')+1:linkFiles[0].find('-')]

for faiFile in faiFiles:
    if querySpecies in faiFile:
        mainSpeciesFai = faiFile

listKaryotypes = fai2karyotype(faiFilePath,karyotypesFilesPath,faiFiles,mainSpeciesFai)


for linkFile in linkFiles:
    KaryotypeFiles = []
    speciesList = [0,0]
    speciesList[0] = linkFile[linkFile.find('.')+1:linkFile.find('-')]
    speciesList[1] = linkFile[linkFile.find('-')+1:linkFile.rfind('_')].replace('PAC2_0.','').replace('PAC4GC.','')
    for species in speciesList:
        for karyotype in listKaryotypes:
            if species in karyotype:
                KaryotypeFiles.append(karyotype)
    if '523' in speciesList and '524' in speciesList:
        subgenomesCheck = 1
    else:
        subgenomesCheck = 0
    generateConfigs(circosConfigFilesPath,karyotypesFilesPath,LinkPath,KaryotypeFiles,linkFile,subgenomesCheck)
    subprocess.call(['/Applications/circos-0.69/bin/circos','-conf',circosConfigFilesPath+'circos.conf','-outputfile',
                     '%s-%s'%(speciesList[0],speciesList[1]),'-outputdir',circosOutPath])
    #print '/Applications/circos-0.69/bin/circos -conf %s'%(configPath+'circos.conf')
    #