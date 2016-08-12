import statistics
speciesHistogramsDict = {'524':[],'523':[],'308':[],'312':[],'313':[],'314':[],'323':[],'316':[]}

inputBed = open('finalSyntenyMultipleSpecies.bed','r')

#524_Chr01K	4799099	5265487	308_Chr01_2552442_2842086|312_scaffold_1_7299767_7652320

for line in inputBed:
    lineList = line.split()
    lineList2 = lineList[-1].split('|')

    speciesHistogramsDict[lineList[0].split('_')[0]].append(int(lineList[2])-int(lineList[1]))

    for item in lineList2:
        if int(item.split('_')[-1])-int(item.split('_')[-2]) >0:
            speciesHistogramsDict[item.split('_')[0]].append(int(item.split('_')[-1])-int(item.split('_')[-2]))

for item in speciesHistogramsDict:
    print item, statistics.mean(speciesHistogramsDict[item])

