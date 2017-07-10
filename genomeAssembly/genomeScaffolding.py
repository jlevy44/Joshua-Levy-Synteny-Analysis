import subprocess, os, sys
from collections import defaultdict, OrderedDict
import numpy as np
from multiprocessing import Pool, Queue, Process
from threading import Thread
import subprocess,shutil
from pybedtools import BedTool
from jcvi.formats import gff
from pyfaidx import Fasta
import time


"""python genomeScaffolding.py ReferenceBuild sampleBuild  CDSProtID OldCDSGeneName protID1 weight1 protID2 weight2 ..."""

CDSgeneNaming = sys.argv[4]
CDSspecies = sys.argv[3]
args = sys.argv[5:]
root = os.getcwd()+'/'
weights = OrderedDict()
listSamplesv0 = [folder for folder in os.listdir('v0') if folder.endswith('v0')]

try:
    ReferenceBuild = int(sys.argv[1])
except:
    ReferenceBuild = 1

try:
    sampleBuild = int(sys.argv[2])
except:
    sampleBuild = 1




print args
print CDSgeneNaming
print CDSspecies

for i in np.arange(0,len(args),2):
    try:
        weights[args[i]]=int(args[i+1])
    except:
        print args
print weights

runCommand = lambda x: subprocess.call(x,shell=True)

binbash = "#!/bin/bash"
makeTrashFolder = 'mkdir oldFiles'
moduleLoads = """module load cufflinks/2.2.1
module load samtools/1.3.1
module load gmap
module load parallel/20150222
module load bedtools/2.25.0
module unload gcc
module load gcc/6.3.0
"""

def runCommands(q):
    while not q.empty():
        print q
        try:
            print q.get()
            runCommand(q.get())
        except:
            with open('Error.txt','a') as f:
                f.write(q.get()+'\n')
        q.task_done()

def buildReferences(reference): # essentially keys of weights
    global root
    global binbash, makeTrashFolder, moduleLoads
    print reference
    os.chdir('./referenceGenomes/'+reference)
    #print os.getcwd()
    #print os.listdir('.')
    fastaOld = [fasta for fasta in os.listdir('.') if 'cds' not in fasta.lower() and (fasta.endswith('.fa') or fasta.endswith('.fasta'))][0]
    #Fasta(fastaOld)
    #gff.load([file for file in os.listdir('.') if 'cufflinks' not in file and (file.endswith('.gff3') or file.endswith('.gff'))][0])
    writeCommands = [binbash,moduleLoads,makeTrashFolder,'samtools faidx %s'%fastaOld,
                     'python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS -o %s'%([file for file in os.listdir('.') if 'cufflinks' not in file and (file.endswith('.gff3') or file.endswith('.gff'))][0],fastaOld,reference+'.cds'),
                     'python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s'%([file for file in os.listdir('.') if 'cufflinks' not in file and (file.endswith('.gff3') or file.endswith('.gff'))][0],reference+'.bed'),
                     'python %sreplacepath.py %s'%(root,reference+'.bed'),'mv %s %s ..'%(reference+'.bed',reference+'.cds')]
    #binbash,makeTrashFolder,moduleLoads,
    #print '\n'.join(writeCommands)
    """if __name__ == '__main__':
        q =  Queue(maxsize=0)
        for command in writeCommands:
            q.put(command)
        runCommands(q)"""
    """for command in writeCommands:
        print command
        try:
            runCommand(command)
        except:
            with open('Error.txt','a') as f:
                f.write(command+'\n')"""
    """for i, command in writeCommands:
        print command
        if (i == 3 or i==4) and (reference + '.bed' not in os.listdir('..') or os.stat('../'+reference + '.bed').st_size == 0):
            runCommand(command)
        elif i == 2 and (reference + '.cds' not in os.listdir('..') or os.stat('../'+reference + '.cds').st_size == 0):
            runCommand(command)
        elif i not in range(2, 7):
            runCommand(command)"""

    with open('buildReference.sh','w') as f:
        f.write('\n'.join(writeCommands))
    subprocess.call(['nohup','sh','buildReference.sh'])
    os.chdir(root)

#print ReferenceBuild


CDSOld = [fasta for fasta in os.listdir('./referenceGenomes/%s'%CDSspecies) if 'cds' in fasta.lower() and (fasta.endswith('.fa') or fasta.endswith('.fasta'))][0]

linkReferences = ['ln -s %s%s/%s.cds %s.cds\nln -s %s%s/%s.bed %s.bed'%(root,'referenceGenomes',ref,ref,root,'referenceGenomes',ref,ref) for ref in weights.keys()]



def buildSamplesv0(sample): #sample = Bdist_xxx_v0.fa
    global root
    global CDSspecies, CDSOld
    global binbash, makeTrashFolder, moduleLoads
    global CDSgeneNaming, linkReferences
    print sample
    os.chdir('v0/'+sample)
    fastaNew = sample+'.fa'
    geneNaming = sample.replace('_','') # -t is number of worker threads
    runCommand('rm finishBuild.txt')
    writeCommands = [binbash,moduleLoads,makeTrashFolder,'rm -r %s %s.gff3.db %s.chromosome *.iit %s.coords'%(geneNaming,geneNaming,geneNaming,geneNaming),
                     'samtools faidx %s' %fastaNew,
                     'gmap_build --dir=. -d %s %s' % (geneNaming,fastaNew),
                     'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 6 %s > %s 2> %s' % (
                     geneNaming, '../../referenceGenomes/%s/'%CDSspecies + CDSOld, geneNaming + '.gff3', geneNaming + '.log'),
                     'python %srenameGenes.py %s %s %s' %(root,geneNaming + '.gff3', CDSgeneNaming ,geneNaming),
                     'python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s' % (geneNaming + '.gff3', sample + '.bed'),
                     'python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS -o %s' % (
                     geneNaming+'.gff3', fastaNew,sample + '.cds')]+linkReferences+['> finishBuild.txt']

                     #"""'python %sgff2CDSBed.py %s'%(root,geneNaming + '.gff3'),'sortBed -i %s.CDS.bed > %s.CDS2.bed'%(geneNaming,geneNaming),
                     #'python %sformatBed.py s %s v0 1'%(root,geneNaming+'.CDS2'),'bedtools getfasta -name -fi ./%s -bed %s.CDS2.bed -fo %s.cds'%(fastaNew,geneNaming,sample)
                     #]"""#'mv %s %s ..'%(sample+'.cds',sample+'.bed') binbash, moduleLoads, makeTrashFolder,
    #'python -m jcvi.formats.gff load %s %s --feature=CDS --id_attribute=Name -o %s' % (geneNaming + '.gff3', fastaNew,sample + '.cds'),
    #'mergeBed -c 4  -i %s.CDS2.bed > %s.CDS.bed'%(geneNaming,geneNaming)
    #print writeCommands
    #print os.getcwd()
    #open('buildSample.sh', 'w').close()
    """if __name__ == '__main__':
        q = Queue(maxsize=0)
        for command in writeCommands:
            q.put(command)
        runCommands(q)"""
    i=0
    """
    for command in writeCommands:
        #print i,command
        #print i
        if (i == 2 or i == 3 or i == 4) and (geneNaming + '.gff3' not in os.listdir('.') or os.stat(geneNaming + '.gff3').st_size ==0):
            print(command)
            runCommand(command)
        elif i==5 and (sample + '.bed' not in os.listdir('.') or os.stat(sample + '.bed').st_size ==0):
            print(command)
            runCommand(command)
        elif i == 6 and (sample + '.cds' not in os.listdir('.') or os.stat(sample + '.cds').st_size ==0):
            print(command)
            runCommand(command)
        elif i not in range(2,7):
            print(command)
            runCommand(command)
        i+=1
    """
    with open('buildSample.sh', 'w') as f:
        f.write('\n'.join(writeCommands))
    #subprocess.call(['nohup', 'sh', 'buildSample.sh'])
    runCommand('qsub -P plant-analysis.p -N %s -cwd -l high.c -pe pe_slots 16 -e %s %s' % (
    'build'+sample.split('_')[1], 'ErrFile.txt', 'buildSample.sh'))
    while True:
        if os.path.isfile('finishBuild.txt'):
            break
        else:
            time.sleep(10)

    os.chdir(root)

    """try:
        runCommand(command)
    except:
        with open('Error.txt','a') as f:
            f.write(command+'\n')"""
    """with open('buildSample.sh','w') as f:
        f.write('\n'.join(writeCommands))
    try:
        subprocess.call(['nohup','sh','buildSample.sh'])
    except:
        with open('output.txt', 'a') as f:
            f.write('Error in %s'%sample)"""
    """writeCommands2 = [binbash, moduleLoads,'gmap_build --dir=. -d %s %s' % (geneNaming,fastaNew),
                     'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 8 %s > %s 2> %s' % (
                     geneNaming, CDSOld, geneNaming + '.gff3', geneNaming + '.log'),
                          'python %srenameGenes.py %s %s %s' % (root, geneNaming + '.gff3', CDSgeneNaming, geneNaming),
                          'python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s' % (
                          geneNaming + '.gff3', sample + '.bed'),
                          'python -m jcvi.formats.gff bed --type=CDS --key=Name %s -o %s' % (
                              geneNaming + '.gff3', sample + '.CDS.bed'),
                          'bedtools getfasta -name -fi ./%s -bed %s.CDS.bed -fo %s.cds' % (
                          fastaNew, sample, sample)]
        with open('buildSample.sh', 'w') as f:
            f.write('\n'.join(writeCommands2))
        subprocess.call(['nohup', 'sh', 'buildSample.sh'])"""






try:
    os.mkdir('v1')
    for folder in listSamplesv0:
        os.mkdir('v1/%s'%folder.replace('v0','v1'))
        os.mkdir('v1/%s/OldFiles'%folder.replace('v0','v1'))
except:
    pass
buildCorrespondence = {folder:folder.replace('v0','v1') for folder in listSamplesv0}
listSamplesv1 = buildCorrespondence.values()

print listSamplesv1

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

fastaNucOld = [fasta for fasta in os.listdir('./referenceGenomes/%s'%CDSspecies) if 'cds' not in fasta.lower() and (fasta.endswith('.fa') or fasta.endswith('.fasta'))][0]


def generatev1(sample):
    os.chdir('v0/%s'%sample)
    print sample.replace('v0', 'v1')
    global binbash, makeTrashFolder, moduleLoads, root, weights, fastaNucOld, CDSspecies
    #print weights
    print '\n'.join('%s %d'%(key,weights[key]) for key in weights.keys())#weights.keys()#'\n'.join('%s %d'%(key,weights[key]) for key in sorted(weights, key=weights.get, reverse=True).keys())
    print 'hi'
    """if __name__ == '__main__':
        p = Pool(None)
        p.imap(pairwise, [(sample,ref) for ref in weights.keys()])"""
    with open('weights.txt','w') as f:
        f.write('\n'.join([weights.keys()[0]+' %d'%weights[weights.keys()[0]],'%snuc %d'%(CDSspecies,weights[CDSspecies]-1)]+['%s %d'%(key,weights[key]) for key in weights.keys()[1:]]))
    nucCommands = [binbash,moduleLoads]+ ['nucmer -t 6 -p %s %s %s'%(CDSspecies+'nuc',root+'referenceGenomes/%s/'%CDSspecies+fastaNucOld,sample+'.fa'),
                 'delta-filter -m -q -i 85 -u 50 %snuc.delta > %snuc2.delta'%(CDSspecies,CDSspecies),'show-tiling -a %snuc2.delta > %snuc.tiling'%(CDSspecies,CDSspecies)]
    commands1 = [binbash, moduleLoads]+['rm *.anchors *.last *.filtered *.prj']+\
                ['nohup python -m jcvi.compara.catalog ortholog %s %s\nmv %s %s'%(ref,sample,'%s.%s.lifted.anchors'%(ref,sample),'%s.%s.lifted.anchors'%(sample,ref)) for ref in weights.keys()]
    commands2=[binbash, moduleLoads]+['rm multipleMapping.bed','\n'.join('python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --qbed=%s --sbed=%s -o %s'%('%s.%s.lifted.anchors'%(sample,ref),ref+'_syn'+'.bed',sample+'_%ssyn'%ref+'.bed','%s.synteny.bed'%(ref)) for ref in weights.keys()),
                                      'python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --qbed=%s --sbed=%s -o %snuc.synteny.bed'%('nucMap.bed',CDSspecies+'_nucSyn.bed',sample+'_nucSyn.bed',CDSspecies),
         'nohup python -m jcvi.assembly.allmaps mergebed %s -o %s'%(' '.join(['%s.synteny.bed'%(ref) for ref in (weights.keys() + [CDSspecies+'nuc'])]),'multipleMapping.bed')]
    qsub=[binbash,moduleLoads]+['python -m jcvi.assembly.allmaps path --skipconcorde --cpus=32 --ngen=300 --npop=50 multipleMapping.bed %s.fa' % (sample),
         'mv multipleMapping.fasta %sv1/%s/%s.fa' % (root,sample.replace('v0', 'v1'), sample.replace('v0', 'v1'))]
    #'nohup liftOver -gff %s.gff3 multipleMapping.chain %s.gff3 unmapped' % (sample.replace('_',''), sample.replace('_','').replace('v0', 'v1')), ,'mv %s.gff3 ../../v1/%s' % (sample.replace('_','').replace('v0', 'v1'), sample.replace('v0', 'v1'))
    #for ref in weights.keys():
    #    pairwise((sample,ref))
    """if __name__ == '__main__':
        q = Queue(maxsize=0)
        for command in commands:
            q.put(command)
        runCommands(q)"""
    #print '\n'.join(commands)
    with open('nucCommand.sh','w') as f:
        f.write('\n'.join(nucCommands))
    with open('constructv1_1.sh','w') as f:
        f.write('\n'.join(commands1))
    with open('constructv1_2.sh','w') as f:
        f.write('\n'.join(commands2))
    with open('qsub_buildv1.sh','w') as f:
        f.write('\n'.join(qsub))
    print os.listdir('%s/v1/%s'%(root,sample.replace('v0','v1')))
    if '%snuc.tiling'%CDSspecies not in os.listdir('.'):
        runCommand('sh nucCommand.sh')
    #print ['%s.%s.lifted.anchors' %(sample, ref) in os.listdir('.') and os.stat('%s.%s.lifted.anchors' %(sample, ref)).st_size > 0 for ref in weights.keys()]
    print all(['%s.%s.lifted.anchors' %(sample, ref) in os.listdir('.') and os.stat('%s.%s.lifted.anchors' %(sample, ref)).st_size > 0 for ref in weights.keys()]) == 0

    #exit()
    if all([os.path.isfile('%s.%s.lifted.anchors' %(sample, ref)) and os.stat('%s.%s.lifted.anchors' %(sample, ref)).st_size > 0 for ref in weights.keys()]) == 0:
        print sample, ['%s.%s.lifted.anchors' %(sample, ref) in os.listdir('.') and os.stat('%s.%s.lifted.anchors' %(sample, ref)).st_size > 0 for ref in weights.keys()]
        runCommand('sh constructv1_1.sh')
        sampleCount = 0
        for ref in weights.keys():
            sampleCount = replaceGeneNames(sample, ref, sampleCount)
            print 'hello ' + sample, ref

    print 'construct_1' + sample + ' done'

    try:
        tiling2bed('%snuc.tiling'%CDSspecies, CDSspecies, sample, sample+'_%ssyn'%CDSspecies+'.bed')
    except:
        print sys.exc_info()[0]
    #exit()
    print 'hi2'
    replaceGeneNames(sample,CDSspecies,0,1)
    if os.stat('nucMap.bed').st_size == 0:
        exit()
    print 'hi3'
    runCommand('sh constructv1_2.sh')
    try:
        if os.stat('./multipleMapping.bed').st_size > 0:
            runCommand('qsub -P plant-analysis.p -N %s -cwd -l h_rt=50:00:00 -pe pe_slots 32 -e %s %s'%(sample,'ErrFile.txt','qsub_buildv1.sh')) #FIXME pe_slots 16, time limit pe_8
        else:
            with open('ErrFile.txt','a') as f:
                f.write('Multiple Mapping Size 0, unable to build v1...')
    except:
        with open('ErrFile.txt', 'a') as f:
            f.write('Multiple Mapping File does not exist, unable to build v1...')
    os.chdir(root)
    #for command in commands:
    #    print command
    #    runCommand(command)
    #FIXME ADD qsub


def formatSamplev0(sample):
    global root
    commands = ['python %sformatBed.py s %s v0'%(root,sample),'python %sformatCDS.py s %s v0'%(root,sample)]
    for command in commands:
        runCommand(command)
        os.chdir(root)

def formatRef(reference):
    global root
    commands = ['python %sformatBed.py r %s v0' % (root, reference), 'python %sformatCDS.py r %s v0' % (root, reference)]
    for command in commands:
        runCommand(command)
        os.chdir(root)

sampleDist = [listSamplesv0[x:x+7] for x in xrange(0,len(listSamplesv0),7)]
print sampleDist

def buildSampv0List(samplist):
    for sample in samplist:
        try:
            buildSamplesv0(sample)
        except:
            print 'Error building ' + sample

def formatv0List(samplist):
    for sample in samplist:
        try:
            formatSamplev0(sample)
        except:
            print 'Error formatting ' + sample

if __name__ == '__main__':
    with open('output.txt', 'a') as f:
        f.write('Outv1')
    listSamplesv0 = [sample for sample in listSamplesv0 if sample.replace('v0', 'v1') + '.fa' not in os.listdir(
        '%sv1/%s' % (root, sample.replace('v0', 'v1')))]
    print len(listSamplesv0) // 6 + 1
    sampleDist = [listSamplesv0[x:x + len(listSamplesv0) // 6 + 1] for x in
                  xrange(0, len(listSamplesv0), len(listSamplesv0) // 6 + 1)]
    print listSamplesv0
    print sampleDist
    if ReferenceBuild:
        p = Pool(processes=6)
        p.map(buildReferences, weights.keys())
        p.map(func=formatRef, iterable=weights.keys())
        p.close()
        p.join()
    if sampleBuild:
        p = Pool(processes=6)#processes=8
        p.map_async(func=buildSampv0List, iterable=sampleDist)
        p.map_async(func=formatv0List, iterable=sampleDist)
        p.close()
        p.join()
    #for samplelist in sampleDist:
    #    p.map(generatev1, samplelist)
    #for ref in weights.keys():
    #    formatRef(ref)
    #buildReferences('460')
    #formatRef('460')
def reader(q):
    while True:
        sample = q.get()
        try:
            generatev1(sample)

        except:
            print 'Generation Error in ' + sample
            with open('Error.txt', 'a') as f:
                f.write('Generation Error in ' + sample + '\n')

        q.task_done()

def genv1List(samplelist):
    for sample in samplelist:
        #generatev1(sample)
        try:
            generatev1(sample)
        except:
            print 'Error gen v1 in ' + sample

if __name__ == '__main__':
    #for samplelist in sampleDist:
    #q = Queue(maxsize=0)
    #num_threads = 6
    #for i in range(num_threads):
    #    worker = Process(target = reader,args=(q,))
    #    worker.daemon=True
    #    worker.start()
    listSamplesv0 = [sample for sample in listSamplesv0 if sample.replace('v0','v1') + '.fa' not in os.listdir('%sv1/%s'%(root,sample.replace('v0','v1')))]
    print len(listSamplesv0)//6 + 1
    sampleDist = [listSamplesv0[x:x + len(listSamplesv0)//6 + 1] for x in xrange(0, len(listSamplesv0), len(listSamplesv0)//6 + 1)]

    p = Pool()
    p.map_async(genv1List,sampleDist)
    #for sample in samplelist:
    #    p.map(generatev1,args=(sample,))

    p.close()
    p.join()
    #for sample in samplelist:
    #    q.put(sample)
    #q.join()
    """try:
            generatev1(sample)
            break
        except:
            print 'Generation Error in '+ sample
            with open('Error.txt','a') as f:
                f.write('Generation Error in '+ sample + '\n')
                break
                """


"""'gffread -E %s -o- > %s' % (geneNaming + '.gff3', sample + '.cufflinks.gff'),
                          'python %sgff2CDSBed.py %s.cufflinks.gff' % (root, sample),
'gffread -E %s -o- > %s' % (geneNaming + '.gff3', sample + '.cufflinks.gff'),
'gffread -E %s -o- > %s'%([file for file in os.listdir('.') if 'cufflinks' not in file and (file.endswith('.gff3') or file.endswith('.gff'))][0],reference+'.cufflinks.gff'),
                     'gffread %s -x %s -g %s'%(reference+'.cufflinks.gff',reference+'.cds',fastaOld),
                    'python %sgff2CDSBed.py %s.cufflinks.gff'%(root,sample),
                     'bedtools getfasta -name -fi ./%s -bed %s.cufflinks.CDS.bed -fo %s.cds'%(fastaNew,sample,sample), """