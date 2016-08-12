pathOfComparison = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/' \
                   'DNA_Segment_Extraction/OldSortFiles/'
open('ComparisonOutput.txt','w').close()
createListFile = open('SortComparisonList.txt','r')


sortFileList=[]

for line in createListFile:
    sortFileList+=[line.split()]

createListFile.close()

comparisonOutput = open('ComparisonOutput.txt','w')

# comparing new sort files to old sort files for all sort files listed
for sortFile in sortFileList:
    speciesNumber = sortFile[1]
    newSortFile = open(sortFile[0],'r')
    oldSortFile = open(pathOfComparison+sortFile[0],'r')
    # create lists of genes in new and old sort files
    newGeneList=[]
    oldGeneList=[]
    comparisonOutput.write('\n\nHere is the comparison between new and old sort file for species '+speciesNumber+':\n')
    for line in newSortFile:
        newGeneList += [line.split()[0]]
    for line in oldSortFile:
        oldGeneList += [line.split()[0]]
    comparisonOutput.write('The following are genes that are in the new sort but not in old sort:\n')
    for entry in newGeneList:
        if entry not in oldGeneList:
            comparisonOutput.write(entry+'\n')
    comparisonOutput.write('The following are genes that are in the old sort but not in new sort:\n')
    for entry in oldGeneList:
        if entry not in newGeneList:
            comparisonOutput.write(entry+'\n')
    newSortFile.close()
    oldSortFile.close()

comparisonOutput.close()
