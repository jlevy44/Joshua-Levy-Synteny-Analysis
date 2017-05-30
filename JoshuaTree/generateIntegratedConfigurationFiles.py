def generateIntAnalysisConfig(analysisName,outputTuple):
    """generateIntAnalysisConfig inputs name of analysis to perform and relevant information in output tuple and
    generates configuration files for downstream analyses."""

    # if performing multiple synteny analysis using SyntenyFinal
    if analysisName == 'syntenyAnalysis':
        open('syntenicConfig.txt','w').close()
        syntenicConfigFile = open('syntenicConfig.txt','w')

        syntenicConfigText = """NameAnalysis = %s
writeFastaOut = %s
Loci_Threshold = %s
pathPythonModules = %s
pathUnOut = %s
pathSort = %s
BPsMergeDist = %s
softMasked = %s
%s Test Analysis path = %s
%sStop


Fasta
genomePath = %s
listOfGenomeFiles:
%sStop"""%outputTuple
        syntenicConfigFile.write(syntenicConfigText)
        syntenicConfigFile.close()

    # if wishing to generate circos graphs
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

    # if choosing to reconstruct a genome using allmaps and large syntenic segments, might not be as reliable as the
    #next option specified
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

    # does the same thing as the previous option, but instead of using large syntenic segments, uses small syntenic
    # genes to generate more data points for accurate build, recommended that this option is selected
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