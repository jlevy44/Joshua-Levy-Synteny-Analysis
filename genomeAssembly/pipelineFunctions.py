from pybedtools import BedTool
from collections import defaultdict, OrderedDict
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


def replaceGeneNames(sample,ref,count=0,nuc=0):
    refGeneCount = 0
    synmap = '%s.%s.lifted.anchors' % (sample, ref)
    if nuc:
        nucAdd = 'nuc'
        synmap = 'nucMap.bed'
        refbed = ref + '_nucSyn.bed'
        sampbed = sample + '_nucSyn.bed'
        a, b = 1, 0
    else:
        nucAdd = ''
        refbed = ref + '.bed'
        sampbed = sample + '.bed'
        a, b = 0, 1
    sampleProt = sample.split('_')[1]
    with open(refbed,'r') as f:
        refBedLines = f.readlines()
    refBedOut = []
    refGenes = defaultdict(list)
    for line in refBedLines:
        if line:
            refGenes[line.split('\t')[3]] = ref+nucAdd+'_'+str(refGeneCount)
            refBedOut.append(line.replace(line.split('\t')[3],ref+nucAdd+'_'+str(refGeneCount)))
            refGeneCount+=1
    #ref+'_syn'+'.bed',sample+'_%ssyn'%ref+'.bed'
    #print refGenes
    with open(sampbed,'r') as f:
        sampBedLines = f.readlines()
    sampBedOut = []
    sampGenes = defaultdict(list)
    for line in sampBedLines:
        if line:
            sampGenes[line.split('\t')[3]] = sampleProt+nucAdd+'_'+str(count)
            sampBedOut.append(line.replace(line.split('\t')[3], sampleProt + nucAdd + '_' + str(count)))
            count+=1
    with open(synmap,'r') as f:
        synRead = f.readlines()
    synOut = []
    for line in synRead:
        if line and '###' not in line:
            try:
                genes = line.split('\t')
                print genes
                synOut.append(line.replace(genes[0],refGenes[genes[a]]).replace(genes[1],sampGenes[genes[b]]))
            except:
                with open('Err.txt','a') as f:
                    f.write(line+'\n')
    """
    if nuc:
        print sampBedOut[0:10]
        print refBedOut[0:10]
        print sampGenes.items()[0:10]
        print refGenes.items()[0:10]
        print synOut[0:10]
    with open('nucMap.bed','r') as f:
        print f.readlines()[0:10]
    """
    if nuc == 0:
        for writeTuple in [(ref+'_syn'+'.bed',refBedOut),(sample+'_%ssyn'%ref+'.bed',sampBedOut),(synmap,synOut)]:
            with open(writeTuple[0],'w') as f:
                f.writelines(writeTuple[1])
    else:
        for writeTuple in [(refbed,refBedOut),(sampbed,sampBedOut),(synmap,synOut)]:
            with open(writeTuple[0],'w') as f:
                f.writelines(writeTuple[1])
    return count

def tiling2bed(tilingFile,ref,sample,sampBed):
    with open(tilingFile,'r') as f:
        tilingLines = f.read().split('\n')
    genesDict = defaultdict(list)
    with open(ref+'_nucSyn.bed','w') as f1, open(sample+'_nucSyn.bed','w') as f2:
        for line in tilingLines:
            if line:
                lineList = line.split('\t')
                int1 = sorted(map(int,lineList[0:2]))
                int1[0] -= 1
                int2 = sorted(map(int,lineList[2:4]))
                int2[0] -= 1
                f1.write('\t'.join([lineList[-2]]+map(str,int1)+['_'.join([lineList[-2]]+map(str,int1)),'0','+']) + '\n')
                f2.write('\t'.join([lineList[-1]]+map(str,int2)+['_'.join([lineList[-1]]+map(str,int2)),'0','+']) + '\n')
                genesDict['_'.join([lineList[-1]]+map(str,int2))] = '_'.join([lineList[-2]]+map(str,int1))
    b = BedTool(sample+'_nucSyn.bed').subtract(BedTool(sampBed),A=True)
    #print b.head()
    #print genesDict.keys()[0:10]
    origGenes = set(genesDict.keys())
    #print str(b).split('\n')[0:10]
    #print [ line.split('\t')[3] for line in str(b).split('\n') if line][0:10]
    remainGenes = set([ line.split('\t')[3] for line in str(b).split('\n') if line])
    #print list(remainGenes)[0:10]
    BadGenes = list(origGenes - remainGenes)
    #print BadGenes[0:10]
    #print len(origGenes), len(remainGenes), len(BadGenes)
    #exit()
    for gene in BadGenes:
        try:
            del genesDict[gene]
        except:
            pass
    with open('nucMap.bed','w') as f:
        f.write('\n'.join('%s\t%s\t100'%item for item in genesDict.items() if item))