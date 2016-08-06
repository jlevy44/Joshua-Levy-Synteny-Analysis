def generateIntAnalysisConfig(analysisName,outputTuple):


    if analysisName == 'syntenyAnalysis':
        open('syntenicConfig.txt','w').close()
        syntenicConfigFile = open('syntenicConfig.txt','w')

        syntenicConfigText = """NameAnalysis = %s
writeFastaOut = %s
Loci_Threshold = %s
pathPythonModules = %s
pathUnOut = %s
pathSort = %s
%s Test Analysis path = %s
%sStop


Fasta
genomePath = %s
listOfGenomeFiles:
%sStop"""%outputTuple
        syntenicConfigFile.write(syntenicConfigText)
        syntenicConfigFile.close()

    if analysisName == 'circosAnalysis':
        open('circosFiguresPipelineConfig.txt', 'w').close()
        circosConfigFile = open('circosFiguresPipelineConfig.txt','w')

        circosConfigText = """# Circos Figures Config File Joshua Levy
faiFiles
%sStop
BedFiles
%sStop
faiFilePath = %s
karyotypesFilesPath = %s
BedInputPath = %s
circosConfigFilesPath = %s
LinkPath = %s
circosOutPath = %s
pathPythonModules = %s
BPsThreshold = %s"""%outputTuple
        circosConfigFile.write(circosConfigText)
        circosConfigFile.close()

    if analysisName == 'reconstructGenome':
        open('%sallMAPConfig.txt'%outputTuple[0],'w').close()
        allmapConfigFile = open('%sallMAPConfig.txt'%outputTuple[0],'w')
        allmapConfigText = """# Allmaps configuration file
pythonPath = %s
systemPath = %s
bedPath = %s
allMAPImageOutputPath = %s
fastaInputName = %s
fastaOutputName = %s
BedFileList
%sStop"""%outputTuple[1:]
        allmapConfigFile.write(allmapConfigText)
        allmapConfigFile.close()

    if analysisName == 'reconstructGenome2':
        open('%sallMAPConfig2.txt' % outputTuple[0], 'w').close()
        allmapConfigFile = open('%sallMAPConfig2.txt' % outputTuple[0], 'w')
        allmapConfigText = """# Allmaps configuration file
pythonPath = %s
systemPath = %s
bedPath = %s
allMAPImageOutputPath = %s
fastaInputName = %s
fastaOutputName = %s
UnOutPath = %s
SortPath = %s
loci_threshold = %s"""%outputTuple[1:]
        allmapConfigFile.write(allmapConfigText)
        allmapConfigFile.close()