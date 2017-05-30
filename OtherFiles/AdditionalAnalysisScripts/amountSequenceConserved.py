import os
import numpy as np
speciesDict = {'ABR113D':'ABRD','ABR113S':'ABRS','Bhyb36S':'BhybsixS','Bhyb36D':'BhybsixD','Bdistachyon':'Bdistachyon','Bstacei':'Bstacei'}
speciesList = ['ABR113D','ABR113S','Bhyb36S','Bhyb36D']
path = '/global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/CNSAnalysisOutputs/'
path2 = '/ConservedElementsBedFiles/'
analyses = ('AllCNSElements','ConservedElements')
def amountConserved(speciesList,species2,analysisType): # species
    global speciesDict, path, path2
    # species ABR113D  ABR113S  Bhyb36D  Bhyb36S
    # analysisType = AllCNSElements, ConservedElements
    finalDict = {'ABR113D':[],'ABR113S':[],'Bhyb36S':[],'Bhyb36D':[]}
    for species in speciesList:
        filename = path+species+path2+[file for file in os.listdir(path+species+path2) if speciesDict[species2] in file and analysisType in file][0]
        with open(filename,'r') as f:
            finalDict[species] = np.sum(np.vectorize(lambda line: int(line.split('\t')[2])-int(line.split('\t')[1]))([line for line in f if line]))
            f.close()
        #print filename, final
    return finalDict
sequenceConserved = {key:{key3:{key2:0 for key2 in speciesList} for key3 in speciesDict.keys()} for key in analyses}#dict.fromkeys(analyses,dict.fromkeys(tuple(speciesList),dict.fromkeys(tuple(speciesDict.keys()),0)))
for analysis in analyses:
    for species2 in speciesList:
        for species1 in speciesDict.keys():
            with open(path+species2+path2+[filename for filename in os.listdir(path+species2+path2) if speciesDict[species1] in filename and analysis in filename][0],'r') as f:
                sequenceConserved[analysis][species1][species2] = np.sum(np.vectorize(lambda line: int(line.split('\t')[2])-int(line.split('\t')[1]))([line for line in f if line]))
                f.close()
#print sequenceConserved
"""
#with open('output.txt','w') as f2:
print sequenceConserved
print speciesDict.keys()
for analysis in analyses:
    print sequenceConserved
    for species in speciesDict.keys():
        print species
        sequenceConserved[analysis][species] = amountConserved(speciesList,species,analysis)
        print sequenceConserved[analysis][species]
            #f2.write('\t'.join(item for item in (analysis,species,species2,str(sequenceConserved[analysis][species][species2]))))
    #f2.close()
print sequenceConserved
buildDict1 = {key: ['','',0,0,0] for key in speciesDict.keys()}
buildDict2 = {key:buildDict1 for key in speciesList}
sequenceConservedOutput = {key:buildDict2 for key in analyses}#dict.fromkeys(analyses,dict.fromkeys(tuple(speciesDict.keys()),['','',0,0,0]))
for analysis in analyses:
    for species in speciesDict.keys():
        sequenceAmount = np.vectorize(lambda x: sequenceConserved[analysis][species][x])(speciesList)
        sequenceAmount2 = [(species2,str(sequenceConserved[analysis][species2][species]))for species2 in speciesList] #np.vectorize(lambda x: sequenceConserved[analysis][x][species])(speciesList)
        sequenceConservedOutput[analysis][species] = [analysis,species,str(np.sum(sequenceAmount)/4),str(np.std(sequenceAmount)),'OriginalAmounts-'+','.join('%s:%s'%tuple(amount) for amount in sequenceAmount2)]#str(amount) for amount in sequenceAmount)]
finalOutput = open('output.txt','w')
for analysis in analyses:
    for species in speciesDict.keys():
        finalOutput.write('\t'.join(item for item in sequenceConservedOutput[analysis][species])+'\n')"""
with open('output.txt','w') as f:
    for analysis in analyses:
        for species in speciesDict.keys():
            #print str(np.mean(sequenceConserved[analysis][species].items()))
            #print str(np.std(sequenceConserved[analysis][species].items()))
            #print ','.join(key+':'+str(sequenceConserved[analysis][species][key]) for key in sequenceConserved[analysis][species].keys())
            f.write('%s\n\n'%('\t'.join(item for item in ("Analysis="+analysis,"Species="+species,"MeanSequence="+str(np.mean(sequenceConserved[analysis][species].values())),
                    "StD="+str(np.std(sequenceConserved[analysis][species].values())),"Originals="+','.join(key+':'+str(sequenceConserved[analysis][species][key]) for key in sequenceConserved[analysis][species].keys())))))
    f.close()