import os
os.chdir('../../..')
from pipelineFunctions import parseConfigFindList, parseConfigFindPath
root = os.getcwd()+'/'
print root


with open('scaffoldConfig.txt','r') as f:
    (version,buildSample,buildRef,constructSample,CDSspecies,CDSOld,CDSgeneNaming, BB, nuc) = tuple([parseConfigFindPath(x,f) for x in
                                                        ['version','buildSample','buildRef','constructSample','CDS','CDSFasta','geneNameOld','com_bb','com_Nuc']])
    weights = parseConfigFindList('weights',f)

try:
    BB = int(BB)
except:
    BB = 0

try:
    nuc = int(nuc)
except:
    nuc = 0

weightsText = '\n'.join([weights[0]]+[item for item in [CDSspecies+'nuc ' + str(int(weights[1].split()[-1]))] if nuc]+[item for item in [CDSspecies+'BB ' + str(int(weights[1].split()[-1]))] if BB]+weights[1:])
weights = {line.split()[0]:int(line.split()[-1]) for line in weights}

with open('references.txt','w') as f:
    f.write('['+','.join("'%s'"%ref for ref in weights.keys())+']')

listSamples = [folder for folder in os.listdir(version) if folder.endswith(version)]
#listSamples = ['Bdist_100_v0','Bdist_001_v0','Bdist_011_v0']


headSh = """#!/bin/bash
module load cufflinks/2.2.1
module load samtools/1.3.1
module load gmap
module load parallel/20150222
module load bedtools/2.25.0
module unload gcc
module load gcc/6.3.0
module load bbtools
"""

print weightsText

for reference in weights.keys():
    print reference
    os.chdir('./referenceGenomes/' + reference)
    fastaOld = [fasta for fasta in os.listdir('.') if
                'cds' not in fasta.lower() and (fasta.endswith('.fa') or fasta.endswith('.fasta'))][0]
    with open('buildRef.sh','w') as f:
        f.write('\n'.join([headSh,'samtools faidx %s' % fastaOld,
            'python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS -o %s' % (
            [file for file in os.listdir('.') if 'cufflinks' not in file and reference in file and (file.endswith('.gff3') or file.endswith('.gff'))][
                0], fastaOld, reference + '.cds'),
            'python -m jcvi.formats.gff bed --type=mRNA --key=Name %s -o %s' % (
            [file for file in os.listdir('.') if 'cufflinks' not in file and reference in file and (file.endswith('.gff3') or file.endswith('.gff'))][
                0], reference + '.bed'),
            'python %sreplacepath.py %s' % (root, reference + '.bed'), 'mv %s %s ..' % (reference + '.bed', reference + '.cds')]+
                          ['cd '+root,'python %sformatBed.py r %s %s' % (root, reference,version),'cd '+root, 'python %sformatCDS.py r %s %s' % (root, reference,version)]))
    os.chdir(root)

linkReferences = ['ln -s %s%s/%s.cds %s.cds\nln -s %s%s/%s.bed %s.bed'%(root,'referenceGenomes',ref,ref,root,'referenceGenomes',ref,ref) for ref in weights.keys()]
fastaNucOld = [fasta for fasta in os.listdir('./referenceGenomes/%s'%CDSspecies) if 'cds' not in fasta.lower() and (fasta.endswith('.fa') or fasta.endswith('.fasta'))][0]

nextVersion = 'v' + str(int(version.strip('v').strip('\n'))+1)

for sample in listSamples:
    print sample.replace(version, nextVersion)

    os.chdir(version + '/' + sample)
    with open('weights.txt','w') as f:
        f.write(weightsText)
    fastaNew = sample + '.fa'
    geneNaming = sample.replace('_', '')
    writeBuild = [headSh,'rm -r %s %s.gff3.db %s.chromosome *.iit %s.coords' % (geneNaming, geneNaming, geneNaming, geneNaming),
        'samtools faidx %s' % fastaNew,
        'gmap_build --dir=. -d %s %s' % (geneNaming, fastaNew),
        'gmap --dir=. -d %s -B 5 -A --format=gff3_gene -n 1 -t 6 %s > %s 2> %s' % (
            geneNaming, '../../referenceGenomes/%s/' % CDSspecies + CDSOld, geneNaming + '.gff3', geneNaming + '.log'),
        'python %srenameGenes.py %s %s %s' % (root, geneNaming + '.gff3', CDSgeneNaming, geneNaming),
        'python %sfixGFFCoordinates.py %s'%(root,geneNaming),
        'python -m jcvi.formats.gff bed --type=mRNA --key=gene_name %s -o %s' % (geneNaming + '.gff3', sample + '.bed'),
        'gffread -x %s.cds -g %s.fa %s.gff3 -E'%(sample,sample,geneNaming)]+linkReferences + ['cd '+root,'python %sformatBed.py s %s %s'%(root,sample,version),'cd '+root,'python %sformatCDS.py s %s %s'%(root,sample,version)]
    """'python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS -o %s' % (
                 geneNaming+'.gff3', fastaNew,sample + '.cds')""" # key=Name was original, now gene_name with
    nucCommands = [headSh]+ ['nucmer -t 6 -p %s %s %s'%(CDSspecies+'nuc',root+'referenceGenomes/%s/'%CDSspecies+fastaNucOld,sample+'.fa'),
                 'delta-filter -m -q -i 85 -u 50 %snuc.delta > %snuc2.delta'%(CDSspecies,CDSspecies),'show-tiling -a %snuc2.delta > %snuc.tiling'%(CDSspecies,CDSspecies)]
    bbCommands = [headSh.replace('module load samtools/1.3.1\n','module unload samtools\nmodule load samtools/1.4\n')] + ['rm -r ref\nrm BBmapped.bed'] + ['bbmap.sh fastareadlen=600 in=%s.fa ref=%s minid=0.97 ef=0.01 outm=BBmapped.bam ambiguous=toss'%(sample,root+'referenceGenomes/%s/'%CDSspecies+fastaNucOld),
                             'python -m jcvi.formats.sam bed BBmapped.bed BBmapped.bam']#threads=6
    commands1 = [headSh]+['rm *.anchors *.last *.filtered *.prj']+\
                ['nohup python -m jcvi.compara.catalog ortholog %s %s\nmv %s %s'%(ref,sample,'%s.%s.lifted.anchors'%(ref,sample),'%s.%s.lifted.anchors'%(sample,ref)) for ref in weights.keys()]
    commands2=[headSh]+['rm multipleMapping.bed','\n'.join('python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --qbed=%s --sbed=%s -o %s'%('%s.%s.lifted.anchors'%(sample,ref),ref+'_syn'+'.bed',sample+'_%ssyn'%ref+'.bed','%s.synteny.bed'%(ref)) for ref in weights.keys())] \
              + [item for item in ['python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --qbed=%s --sbed=%s -o %snuc.synteny.bed'%('nucMap.bed',CDSspecies+'_nucSyn.bed',sample+'_nucSyn.bed',CDSspecies)] if nuc] \
              + [item for item in ['python -m jcvi.assembly.syntenypath bed %s --switch --scale=10000 --qbed=%s --sbed=%s -o %sBB.synteny.bed' % (
                        'BBMap.bed', CDSspecies + '_BBSyn.bed', sample + '_BBSyn.bed', CDSspecies)] if BB] \
              + ['nohup python -m jcvi.assembly.allmaps mergebed %s -o %s'%(' '.join(['%s.synteny.bed'%(ref) for ref in (weights.keys() + [item for item in [CDSspecies+'nuc'] if nuc] + [item for item in [CDSspecies+'BB'] if BB])]),'multipleMapping.bed')]
    qsub=[headSh]+['python -m jcvi.assembly.allmaps path --skipconcorde --cpus=16 --ngen=400 --npop=60 multipleMapping.bed %s.fa' % (sample),
         'mv multipleMapping.fasta %s%s/%s/%s.fa' % (root,nextVersion,sample.replace(version, nextVersion), sample.replace(version, nextVersion))]

    with open('build.sh','w') as f:
        f.write('\n'.join(writeBuild))
    with open('nucCommand.sh','w') as f:
        f.write('\n'.join(nucCommands))
    with open('constructv1_1.sh','w') as f:
        f.write('\n'.join(commands1))
    with open('constructv1_2.sh','w') as f:
        f.write('\n'.join(commands2))
    with open('BB_build.sh','w') as f:
        f.write('\n'.join(bbCommands))
    with open('qsub_build.sh','w') as f:
        f.write('\n'.join(qsub))

    os.chdir(root)



