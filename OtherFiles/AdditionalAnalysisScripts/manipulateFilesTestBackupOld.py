import shutil, subprocess
import numpy as np
import re
import sys

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname
            configFile.seek(0)
            return line.split()[-1].strip('\n')
configFile = open('configCNSAnalysis.txt','r')
rootfolder = parseConfigFindPath('root_folder',configFile)
pathPython = parseConfigFindPath('pathPython',configFile)
configFile.close()

sys.path.append(pathPython)

outFname = ['NKHDS11111.txt','NKHDS01111.txt','NKHDS10111.txt','NKHDS11100.txt','NKHDS01100.txt','NKHDS10100.txt','NKH111.txt','NKH011.txt','NKH101.txt','HDS111.txt']
paths = ['%s/CS/NKH111/'%rootfolder,
         '%s/CS/NKH011/'%rootfolder,
         '%s/CS/NKH101/'%rootfolder,
         '%s/CS/HDS111/'%rootfolder]

def bed2dict(bedfname):
    bedDict={}
    for line in open(bedfname,'r'):
        if line:
            MASeg = line.split()[-1].split(';')[2].replace('SegmentID=','')
            if bedDict.has_key(MASeg) and tuple(line.split()[0:3]) not in bedDict[MASeg]:
                bedDict[MASeg].append(tuple(line.split()[0:3]))
            else:
                bedDict[MASeg] = [tuple(line.split()[0:3])]
    return bedDict

def grabPartMAFSeq(MAFSeq,x_coords):
    coordtransform = np.array(range(len(MAFSeq)))
    coordtransform = np.delete(coordtransform,[val.start(0) for val in re.finditer('-',MAFSeq)])
    print coordtransform, coordtransform[x_coords[1]-1], x_coords[1]
    return MAFSeq[coordtransform[x_coords[0]]:coordtransform[x_coords[1]-1]]


listALLFiles = str(subprocess.Popen('ls', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                   .stdout.read()).split('\n')

ultimateBedDicts = {'Conserved_CDS.bed':{},'CS_Minus_CDS.bed':{}}

for file in listALLFiles:
    if file.endswith('Conserved_CDS.bed') or file.endswith('CS_Minus_CDS.bed'):
        ultimateBedDicts[file[file.find('_')+1:]][file[:file.find('_')]] = bed2dict(file)

for path in paths:
    for file in str(subprocess.Popen(['ls',path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                              .stdout.read()).split('\n'):
        if file and file.endswith('.fasta'):
            """>Pvirgatum.Pvirgatum.Chr01K_14565140_14836634_328_23_+
            AGACAGCTAACAGCTTACTTTCT
            >Phallii.Phallii.Chr03_5677168_5837565_290_23_+
            AGACAGCTAACAGCTTACTTTCT
            """
            inputFasta=open(path+file,'r') #011
            MASeg = file.split('.')[0]
            speciesList = []
            for line in inputFasta:
                if '>' in line:
                    lineList = line[1:].split('.')
                    speciesName = lineList[0]
                    lineList2 = lineList[-1].split('_')
                    scaffold_name = lineList[-1][:lineList[-1].find(lineList2[-5])-1]
                    #orientation = lineList2[-1].strip('\n')
                    if 'K' == scaffold_name[-1]:
                        speciesName += 'K'
                    if 'N' == scaffold_name[-1]:
                        speciesName += 'N'
                    """
                    if orientation == '+':
                        x = (int(lineList2[-5])+int(lineList2[-3]),int(lineList2[-5])+int(lineList2[-3])+int(lineList2[-2]))
                        print (scaffold_name,),x
                    else:
                        x = (int(lineList2[-4])-int(lineList2[-3])-int(lineList2[-2]),int(lineList2[-4])-int(lineList2[-3]))
                        print (scaffold_name,),x
                    """
                else:
                     speciesList.append(speciesName)
                    #for key in ultimateBedDicts.keys():
                    #    if ultimateBedDicts[key][speciesName].has_key(MASeg):
                    #        for CNSCDS in ultimateBedDicts[key][speciesName][MASeg]:
                                #xRel = tuple(np.vectorize(lambda y: int(y))(CNSCDS[1:])-x[0])
                                #print ultimateBedDicts[key][speciesName][MASeg],conservedSeg, len(conservedSeg), xRel,grabPartMAFSeq(conservedSeg,xRel)
            CDS_conditional = np.vectorize(lambda speciesName: ultimateBedDicts['Conserved_CDS.bed'][speciesName].has_key(MASeg))(speciesList)
            if np.all(CDS_conditional):
                shutil.copy(path+file,path.replace('CS','ProteinCoding_CS')+file)
            elif np.any(CDS_conditional) and not np.all(CDS_conditional):
                shutil.copy(path + file, path.replace('CS', 'Mixed_CS') + file)
            elif not np.any(CDS_conditional):
                shutil.copy(path + file, path.replace('CS', 'NonProteinCoding_CS') + file)
            inputFasta.close()

def concatFastas(fname,HDSPath,NKHPath,FinalPath):
    NKHFileStringList = open(NKHPath+fname,'r').read().split('>')
    HDSFileStringList = open(HDSPath+fname,'r').read().split('>')
    outputStringList = list(set(NKHFileStringList+HDSFileStringList))
    outputFile = open(FinalPath+fname,'w')
    outputFile.write('>'.join(seqs for seqs in outputStringList))
    outputFile.close()

def grabConservedSeq(PathAndfname,conservedSeqs):
    inputFastaReadSegments = open(PathAndfname,'r').read().split('>')
    for sequence in inputFastaReadSegments:
        if sequence:
            speciesName = sequence.split('.')[0]
            if 'K' == sequence.split('_')[0][-1]: #FIXME cahnge _ to -
                speciesName += 'K'
            if 'N' == sequence.split('_')[0][-1]:
                speciesName += 'N'
            if conservedSeqs.has_key(speciesName):
                conservedSeqs[speciesName].append(sequence.split('\n')[1].strip('\n'))
            else:
                conservedSeqs[speciesName] = [sequence.split('\n')[1].strip('\n')]
    return conservedSeqs




for analysisType in ['CS','ProteinCoding_CS','Mixed_CS','NonProteinCoding_CS']:
    arrayFiles = []
    TwoCopyStatNumbers = {}
    analysisPath = paths[0][:-1][:paths[0][:-1].rfind('/')+1].replace('CS',analysisType)
    outFname = ['NKHDS11111.txt','NKHDS01111.txt','NKHDS10111.txt','NKHDS11100.txt','NKHDS01100.txt','NKHDS10100.txt','NKH111.txt','NKH011.txt','NKH101.txt','HDS111.txt']
    for path in paths:
        arrayFiles.append(str(subprocess.Popen(['ls', '%s' % path.replace('CS',analysisType)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        .stdout.read()).split('\n'))
    #print analysisPath
    for i in range(3):
        np.savetxt(analysisPath+outFname[i],np.intersect1d(arrayFiles[3],arrayFiles[i]),fmt='%s',comments='')
        np.savetxt(analysisPath+outFname[i+3],np.setdiff1d(arrayFiles[i],arrayFiles[3]),fmt='%s',comments='')
        np.savetxt(analysisPath+outFname[i+6],arrayFiles[i],fmt='%s',comments='')

    np.savetxt(analysisPath+outFname[9],arrayFiles[3],fmt='%s',comments='')

    folderPaths = {}
    for fname in outFname[0:6]:
        folderPaths[fname.split('.')[0]] = '%s/%s/%s/'%(rootfolder,analysisType,fname.split('.')[0])
        inputFile = open(analysisPath+fname,'r')
        conditionalState = fname.split('.')[0].strip('NKHDS')
        if len(conditionalState) == 5 and conditionalState.endswith('11'):
            HDSpath = '%s/%s/HDS111/'%(rootfolder,analysisType)
            NKHpath = '%s/%s/NKH%s/'%(rootfolder,analysisType,conditionalState[0:3])
            finalPath = '%s/%s/NKHDS%s/'%(rootfolder,analysisType,conditionalState)
            for line in inputFile:
                if line and line != '\n':
                    inputNKHDSfname = line.strip('\n')
                    concatFastas(inputNKHDSfname,HDSpath,NKHpath,finalPath)
        elif len(conditionalState) == 5 and conditionalState.endswith('00'):
            NKHpath = '%s/%s/NKH%s/'%(rootfolder,analysisType,conditionalState[0:3])
            finalPath = '%s/%s/NKHDS%s/'%(rootfolder,analysisType,conditionalState)
            for line in inputFile:
                if line:
                    inputNKHDSfname = line.strip('\n')
                    shutil.copy(NKHpath+inputNKHDSfname,finalPath)
        inputFile.seek(0)
        inputFile.close()
    for fname in outFname:
        fileArray2 = open(analysisPath+fname,'r').read().splitlines()
        speciesRatios = fname.replace('.txt','')
        TwoCopyStatNumbers[speciesRatios] = float(len(fileArray2))
        conservedSeqs = {}
        for file in fileArray2:
            if file:
                conservedSeqs = grabConservedSeq(analysisPath+speciesRatios+'/'+file,conservedSeqs)
        outputFinalFasta = fname.replace('.txt','.fasta')
        open(analysisPath+outputFinalFasta,'w').close()
        outputFinalFastaFile = open(analysisPath+outputFinalFasta, 'w')
        for species in conservedSeqs.keys():
            outputFinalFastaFile.write('>%s\n%s\n'%(species,''.join(conservedSeq for conservedSeq in conservedSeqs[species])))
        outputFinalFastaFile.close()

    ratio2copyStat = [TwoCopyStatNumbers['NKHDS11111']/(TwoCopyStatNumbers['NKHDS10111']+TwoCopyStatNumbers['NKHDS01111']+TwoCopyStatNumbers['NKHDS11111']),
                      TwoCopyStatNumbers['NKHDS11100'] / (TwoCopyStatNumbers['NKHDS10100'] + TwoCopyStatNumbers['NKHDS01100']+TwoCopyStatNumbers['NKHDS11100'])]
    open(analysisPath+'Twocopyratios.txt','w').close()
    twoCopyFile = open(analysisPath+'Twocopyratios.txt','w')
    twoCopyFile.write('Intergenus Conservation [NKHDS11111/(NKHDS11111+NKHDS10111+NKHDS01111)]: %d\n'
                      'Intragenus Conservation [NKHDS11100/(NKHDS11100+NKHDS10100+NKHDS01100)]: %d\n'%(tuple(ratio2copyStat)))
    twoCopyFile.close()