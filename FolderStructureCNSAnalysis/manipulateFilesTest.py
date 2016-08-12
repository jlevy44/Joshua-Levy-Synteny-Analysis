import shutil, subprocess
import numpy as np
import re
import sys, os

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname
            configFile.seek(0)
            return line.split()[-1].strip('\n')
configFile = open('configCNSAnalysis.txt','r')
rootfolder = parseConfigFindPath('root_folder',configFile) # grab root folder for analysis
pathPython = parseConfigFindPath('pathPython',configFile) # grab python path
configFile.close()

sys.path.append(pathPython) # add python modules

# conditionals, will take conditionals in paths, and intersect and differentiate them in different ways to find following combos
# this can find one or two copy (NK) sequences for intra vs inter genus conservation (DS00 or DS11)
outFname = ['NKHDS11111.txt','NKHDS01111.txt','NKHDS10111.txt','NKHDS11100.txt','NKHDS01100.txt','NKHDS10100.txt','NKH111.txt','NKH011.txt','NKH101.txt','HDS111.txt']
paths = ['%s/CS/NKH111/'%rootfolder,
         '%s/CS/NKH011/'%rootfolder,
         '%s/CS/NKH101/'%rootfolder,
         '%s/CS/HDS111/'%rootfolder]

def bed2dict(bedfname): # create a dictionary that has MA Segments as keys
    bedDict={}
    for line in open(bedfname,'r'):
        if line:
            MASeg = line.split()[-1].split(';')[2].replace('SegmentID=','')
            if bedDict.has_key(MASeg) and tuple(line.split()[0:3]) not in bedDict[MASeg]:
                bedDict[MASeg].append(tuple(line.split()[0:3]))
            else:
                bedDict[MASeg] = [tuple(line.split()[0:3])]
    return bedDict

def grabPartMAFSeq(MAFSeq,x_coords): # irrelevant script now...
    coordtransform = np.array(range(len(MAFSeq)))
    coordtransform = np.delete(coordtransform,[val.start(0) for val in re.finditer('-',MAFSeq)])
    print coordtransform, coordtransform[x_coords[1]-1], x_coords[1]
    return MAFSeq[coordtransform[x_coords[0]]:coordtransform[x_coords[1]-1]]


listALLFiles = str(subprocess.Popen('ls', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                   .stdout.read()).split('\n')

ultimateBedDicts = {'Conserved_CDS.bed':{},'CS_Minus_CDS.bed':{}} # split into ~CDS or ~CNS

for file in listALLFiles:
    if file.endswith('Conserved_CDS.bed') or file.endswith('CS_Minus_CDS.bed'):
        ultimateBedDicts[file[file.find('_')+1:]][file[:file.find('_')]] = bed2dict(file) # a dictionary with keys [species] and inside a key of MASEG
        # hierarchy is ultimate[CDS or CNS][species][MASeg]

def concatFastas(fname,HDSPath,NKHPath,FinalPath): # if intergenus conservation, concatenate fastas together and add it to relevant path
    NKHFileStringList = open(NKHPath+fname,'r').read().split('>')
    HDSFileStringList = open(HDSPath+fname,'r').read().split('>')
    outputStringList = list(set(NKHFileStringList+HDSFileStringList))
    outputFile = open(FinalPath+fname,'w')
    outputFile.write('>'.join(seqs for seqs in outputStringList))
    outputFile.close()

def grabConservedSeq(PathAndfname,conservedSeqs): # build a list of conserved sequences that are either Conserved sequence,
    # protein coding, nonprotein, or mixed and will join list together
    inputFastaReadSegments = open(PathAndfname,'r').read().split('>')
    for sequence in inputFastaReadSegments:
        if sequence:
            speciesName = sequence.split('.')[0]
            if 'K' == sequence.split('_')[0][-1]: #FIXME cahnge _ to -
                speciesName += 'K'
            if 'N' == sequence.split('_')[0][-1]:
                speciesName += 'N'
            if conservedSeqs.has_key(speciesName):
                conservedSeqs[speciesName].append(sequence.split('\n')[1].strip('\n')) # add sequence to main list to generate main fasta file output
            else:
                conservedSeqs[speciesName] = [sequence.split('\n')[1].strip('\n')]
    return conservedSeqs




for analysisType in ['CS','ProteinCoding_CS','Mixed_CS','NonProteinCodingN_CS']: # run analysis on these types and sort them into relevant conditionals and output fastas for PhyML analysis and conservation ratios
    # see more details about ratio conservation in various documents...
    TwoCopyStatNumbers = {}
    analysisPath = paths[0][:-1][:paths[0][:-1].rfind('/')+1]#.replace('CS',analysisType)
    outFname = ['NKHDS11111.txt','NKHDS01111.txt','NKHDS10111.txt','NKHDS11100.txt','NKHDS01100.txt','NKHDS10100.txt','NKH111.txt','NKH011.txt','NKH101.txt','HDS111.txt']
    if analysisType == 'CS':
        analysisType2 = ''
        arrayFiles = []
        for path in paths:
            arrayFiles.append(str(subprocess.Popen(['ls', '%s' % path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    .stdout.read()).split('\n'))#.replace('CS',analysisType)
        #print analysisPath
        for i in range(3):
            np.savetxt(analysisPath+outFname[i].replace('.txt','_%s.txt'%analysisType),np.intersect1d(arrayFiles[3],arrayFiles[i]),fmt='%s',comments='')
            np.savetxt(analysisPath+outFname[i+3].replace('.txt','_%s.txt'%analysisType),np.setdiff1d(arrayFiles[i],arrayFiles[3]),fmt='%s',comments='')
            np.savetxt(analysisPath+outFname[i+6].replace('.txt','_%s.txt'%analysisType),arrayFiles[i],fmt='%s',comments='')

        np.savetxt(analysisPath+outFname[9].replace('.txt','_%s.txt'%analysisType),arrayFiles[3],fmt='%s',comments='')
        # the above script finds the intersection and differences
    else:
        analysisType2 = analysisType
        for fname in outFname[0:6]:
            outFile = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'w')
            outFile.write(str(subprocess.Popen('ls %s/%s/%s | grep "%s"'%(rootfolder,'CS',fname.replace('.txt',''),
                    analysisType),shell=True , stdout=subprocess.PIPE).stdout.read()))
            outFile.close() # the CS analysis separates NHK files into various NKHDS folders, now they are named according to type of analysis

    if analysisType == 'CS':
        # if in CS, sort NKH files into various NKHDS folders according to conditional statements and which files fall under each conditional
        for fname in outFname[0:6]:
            inputFile = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'r')
            conditionalState = fname.split('.')[0].strip('NKHDS')
            if len(conditionalState) == 5 and conditionalState.endswith('11'):
                HDSpath = '%s/%s/HDS111/'%(rootfolder,'CS')
                NKHpath = '%s/%s/NKH%s/'%(rootfolder,'CS',conditionalState[0:3])
                finalPath = '%s/%s/NKHDS%s/'%(rootfolder,'CS',conditionalState)
                for line in inputFile:
                    if line and line != '\n' and ':' not in line and '/' not in line:
                        inputNKHDSfname = line.strip('\n')
                        concatFastas(inputNKHDSfname,HDSpath,NKHpath,finalPath)
            elif len(conditionalState) == 5 and conditionalState.endswith('00'):
                NKHpath = '%s/%s/NKH%s/'%(rootfolder,'CS',conditionalState[0:3])
                finalPath = '%s/%s/NKHDS%s/'%(rootfolder,'CS',conditionalState)
                for line in inputFile:
                    if line and line != '\n' and ':' not in line and '/' not in line:
                        inputNKHDSfname = line.strip('\n')
                        #print line, fname.replace('.txt','_%s.txt'%analysisType)
                        #print NKHpath, inputNKHDSfname, finalPath
                        shutil.copy(NKHpath+inputNKHDSfname,finalPath)
            inputFile.seek(0)
            inputFile.close()
    for fname in outFname[0:6]: # compute twocopy statistics (ratio) between NKHDS, see output text
        # output final fasta for NKHDS conditional type and analysis type seen above, which generate trees for each
        fileArray2 = open(analysisPath+fname.replace('.txt','_%s.txt'%analysisType),'r').read().splitlines()
        speciesRatios = fname.replace('.txt','')
        TwoCopyStatNumbers[speciesRatios] = float(len(fileArray2))
        conservedSeqs = {}
        for file in fileArray2:
            if file and ':' not in file and '/' not in file:
                try:
                    conservedSeqs = grabConservedSeq(analysisPath+speciesRatios+'/'+file,conservedSeqs)
                except:
                    print analysisType
        outputFinalFasta = fname.replace('.txt','_%s.fasta'%analysisType)
        open(analysisPath+outputFinalFasta,'w').close()
        outputFinalFastaFile = open(analysisPath+outputFinalFasta, 'w')
        for species in conservedSeqs.keys():
            outputFinalFastaFile.write('>%s\n%s\n'%(species,''.join(conservedSeq for conservedSeq in conservedSeqs[species])))
        outputFinalFastaFile.close()

    ratio2copyStat = [TwoCopyStatNumbers['NKHDS11111']/(TwoCopyStatNumbers['NKHDS10111']+TwoCopyStatNumbers['NKHDS01111']+TwoCopyStatNumbers['NKHDS11111']),
                      TwoCopyStatNumbers['NKHDS11100'] / (TwoCopyStatNumbers['NKHDS10100'] + TwoCopyStatNumbers['NKHDS01100']+TwoCopyStatNumbers['NKHDS11100'])]
    open(analysisPath+'Twocopyratios_%s.txt'%analysisType,'w').close()
    twoCopyFile = open(analysisPath+'Twocopyratios_%s.txt'%analysisType,'w')
    twoCopyFile.write('Intergenus Conservation [NKHDS11111/(NKHDS11111+NKHDS10111+NKHDS01111)]: %f\n'
                      'Intragenus Conservation [NKHDS11100/(NKHDS11100+NKHDS10100+NKHDS01100)]: %f\n'%(tuple(ratio2copyStat)))
    twoCopyFile.close()
    if analysisType == 'CS':
        for path in ['%s/%s/%s/'%(rootfolder,'CS',fname.replace('.txt','')) for fname in outFname[0:6]]:
            for file in str(subprocess.Popen(['ls', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                    .stdout.read()).split('\n'):
                if file and file.endswith('.fasta'):
                    """>Pvirgatum.Pvirgatum.Chr01K_14565140_14836634_328_23_+
                    AGACAGCTAACAGCTTACTTTCT
                    >Phallii.Phallii.Chr03_5677168_5837565_290_23_+
                    AGACAGCTAACAGCTTACTTTCT
                    """
                    inputFasta = open(path + file, 'r')  # 011
                    MASeg = file.split('.')[0]
                    speciesList = []
                    for line in inputFasta:
                        if '>' in line:
                            lineList = line[1:].split('.')
                            speciesName = lineList[0]
                            lineList2 = lineList[-1].split('_')
                            scaffold_name = lineList[-1][:lineList[-1].find(lineList2[-5]) - 1]
                            # orientation = lineList2[-1].strip('\n')
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
                            # for key in ultimateBedDicts.keys():
                            #    if ultimateBedDicts[key][speciesName].has_key(MASeg):
                            #        for CNSCDS in ultimateBedDicts[key][speciesName][MASeg]:
                            # xRel = tuple(np.vectorize(lambda y: int(y))(CNSCDS[1:])-x[0])
                            # print ultimateBedDicts[key][speciesName][MASeg],conservedSeg, len(conservedSeg), xRel,grabPartMAFSeq(conservedSeg,xRel)
                    CDS_conditional = np.vectorize(
                        lambda speciesName: ultimateBedDicts['Conserved_CDS.bed'][speciesName].has_key(MASeg))(speciesList)
                    # in each folder for the NKHDS conditionals, change file names to match relevant analysis, eg. whether
                    # sequence intersects a CDS, is not, or some do and some don't between different species
                    if np.all(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_ProteinCoding_CS.fasta'))
                    elif np.any(CDS_conditional) and not np.all(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_Mixed_CS.fasta'))
                    elif not np.any(CDS_conditional):
                        os.rename(path + file, path + file.replace('.fasta', '_NonProteinCodingN_CS.fasta'))
                    inputFasta.close()