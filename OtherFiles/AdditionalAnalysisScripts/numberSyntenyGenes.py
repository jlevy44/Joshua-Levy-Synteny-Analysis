from pybedtools import *


bedList = open('queryAndComparisonBedList.txt','r')
bedStructure = {}
for line in bedList:
    if 'q' in line:
        bedStructure[line.split('.')[-2]] = [BedTool(line.strip('\n'))]
    else:
        for species in bedStructure.keys():
            if species in line.split('.')[1]:
                bedStructure[species]+=[[BedTool(line.strip('\n')),line.split('.')[-2].split('_')[0]]]

def numberSyntenicGenes(pairedComparison,querySpecies):
    """This function finds number of syntenic genes that exist between paired comparison bed object and species A
    bed object. Output is number..."""
    #NOTE, .sort.merge gets rid of duplicates of overlapping genes
    # getting rid of the above will fail to account for multiple copies of genes (discuss in more detail)
    return pairedComparison.intersect(querySpecies).sort().merge().count()

for species in bedStructure.keys():
    for bedObject in bedStructure[species][1:]:
        print(species,bedObject[1],numberSyntenicGenes(bedObject[0],bedStructure[species][0]))