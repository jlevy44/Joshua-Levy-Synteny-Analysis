import subprocess, os
binbash = "#!/bin/bash"
makeTrashFolder = 'mkdir oldFiles'
moduleLoads = """module load cufflinks/2.2.1
module load samtools/1.3.1
module load gmap
module load parallel/20150222
module load bedtools/2.25.0
module unload gcc
module load gcc/5.4.0
"""
class speciesBuildFiles():
    def __init__(self,fastafile,isReference,referenceSpecies):
        self.fastafile = fastafile
        self.isReference = isReference # should be one if it is a reference species
        self.referenceSpecies = referenceSpecies # species that it will build its [synteny bed, bed, CDS and gff3 file from]
        self.geneNaming = fastafile[:fastafile.find('.')]
        self.cds = self.geneNaming+'.cds'
        self.bed = self.geneNaming+'.bed'

    def buildSample(self,referenceSpecies,idx):
        open('buildSample%d.sh'%idx,'w').close()
        writeCommands = [binbash,moduleLoads,'samtools faidx %s'%fastaOld,'gmap_build --dir=. -d %s %s'%(self.geneNaming,referenceSpecies.cds),
                     'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 8 %s > %s 2> %s'%(self.geneNaming,CDSOld,self.geneNaming+'.gff3',self.geneNaming+'.log'),
                         'python renameGenes.py %s %s' % (self.geneNaming + '.gff3', self.geneNaming),
                             'gffread - E %s - o - > %s'%(self.geneNaming+'.gff3',self.geneNaming+'.cufflinks.gff'),
                     'python gff2CDSBed.py Bhyb.cufflinks.gff','bedtools getfasta -name -fi ./hybridumTest.fa -bed Bhyb.cufflinks.CDS.bed -fo Bhyb.cds',] # FIXME change renamegenes function to take on second argument
        with open('buildSample%d.sh'%idx,'w') as f:
            f.write()




def buildReference(fastaOld,CDSOld): # should be in format Bdist.fa and any other naming scheme for cds file, cds ---> bdist.cds
    speciesName = fastaOld[:fastaOld.find('.')]
    open('buildReference.sh','w').close()
    writeCommands = [binbash,makeTrashFolder,moduleLoads,'samtools faidx %s'%fastaOld,'gmap_build --dir=. -d %s %s'%(speciesName,fastaOld),
                     'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 8 %s > %s 2> %s'%(speciesName,CDSOld,speciesName+'.gff3',speciesName+'.log'),
                             'gffread - E %s - o - > %s'%(speciesName+'.gff3',speciesName+'.cufflinks.gff'),
                     'gffread %s -x %s -g %s'%(speciesName+'.cufflinks.gff',speciesName+'.cds',fastaOld),'python -m jcvi.formats.gff bed %s > %s'%(speciesName+'.gff3',speciesName+'.bed'),
                     'python replacepath.py %s'%speciesName+'.bed']
    with open('buildReference.sh',) as f:
        f.writelines(writeCommands)
    subprocess.call(['nohup','sh','buildReference.sh'])
    return speciesBuildFiles(fastaOld,1,0)






fastaOld =1
CDSOld =1

# doing the below for now
listSamples = [file for file in os.listdir('.') if file and file != fastaOld and file != CDSOld and file.endswith('.fa')] #maybe include .fasta??
reference = buildReference(fastaOld,CDSOld)
sampleBuilds = [speciesBuildFiles(fasta,0,reference) for fasta in listSamples]
[build.buildSample(reference,idx) for idx,build in enumerate(sampleBuilds)]
#throw out unneccesary files into the trash directory