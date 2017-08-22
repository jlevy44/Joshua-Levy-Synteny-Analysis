#!/usr/bin/env nextflow

String findValue(String inputStr) {
configFile = file('kmerCountCompare.txt')
allLines = configFile.readLines()
for( line in allLines ) {
    if(line.startsWith(inputStr)){
        println (line - inputStr - '= ')
        return (line - inputStr - '= ');
    }
}

}

String findPeak(String peakfname){

    String[] parts = peakfname.split("_");
    String peakfinalname = parts[0];
    return peakfinalname;

}

// parse config

pythonPath = findValue('pythonPath ' );
systemPath = findValue('systemPath ' );
blastPath = findValue('blastPath ' );
kmercountPath = findValue('kmercountPath ' );
reclusterPath = findValue('reclusterPath ');
best50kmerPath = findValue('kmer50BestPath ');
fastaPath = findValue('fastaPath ' );
bedPath = findValue('bedPath ' );
sortPath = findValue('sortPath ' );
sortedbedPath = findValue('sortedbedPath ' );
genome = findValue('genome ' );
peakPath = findValue(' peakPath');


//preprocess = findValue('preprocess ').asType(Integer);


save = findValue('save ');
old = findValue('old ');



splitFast = findValue('splitFasta ').asType(Integer);
writeKmer = findValue('writeKmer ').asType(Integer);
fromFasta = findValue('kmer2Fasta ').asType(Integer);
writeBlast = findValue('writeBlast ').asType(Integer);
b2b = findValue('blast2bed ').asType(Integer);
genMat = findValue('generateClusteringMatrix ').asType(Integer);
kHist = findValue('kmerHistogram ').asType(Integer);
bPeak = findValue('blastPeaks ').asType(Integer);
bPeak2 = findValue('peakClusterMatrix ').asType(Integer);
trans = findValue('transformData ').asType(Integer);
trans2 = findValue('transformDataChildren ').asType(Integer);
clust = findValue('ClusterAll ').asType(Integer);



genomeSplitName = genome - '_split' - '.fa' + '_split.fa';
blastDBName = genomeSplitName - '.fa'
genomeFullPath = fastaPath + genomeSplitName;
kmercountName = genomeSplitName - '.fa' + '.kcount' + '.fa';
blastName = kmercountName - '.fa' + '.BLASTtsv.txt';
workingDir = new File('').getAbsolutePath();

check = Channel.from(genomeSplitName,genomeFullPath,kmercountName,blastName,workingDir)
                .subscribe{println it}

//gCh = Channel.from(genome)
//            .subscribe {println it + 'aaaa'}

genomeChan = Channel.from(genome - '\n')

genomeChan2 = Channel.create()

process splitFastaProcess {
executor = 'local'
clusterOptions = { splitFast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan

output:
    stdout genomeChan2

script:
if(splitFast == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py splitFasta ${genomeName} ${fastaPath}
    """
else
    """
    echo ${genomeSplitName}
    """

}

genomeChan3 = genomeChan2.last()

gCh = genomeChan3
                .subscribe{println it -'\n' + 'aaniofs'}

genomeChan4 = Channel.create()

process writeKmerCount {

clusterOptions = { writeKmer == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan3

output:
    val genomeName into genomeChan4


script:
if(writeKmer == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py writeKmerCount ${fastaPath} ${kmercountPath}
    """
else
    """
    touch done
    """

}

genomeChan5 = Channel.create()
genomeChan55 = Channel.create()

process kmer2Fasta {
executor = 'local'
clusterOptions = { fromFasta == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan4

output:
    val genomeName into genomeChan5
    //val genomeName into genomeChan55


script:
if(fromFasta == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    python subgenomeClusteringInterface.py kmer2Fasta ${kmercountPath}
    makeblastdb -in ${genomeFullPath} -dbtype nucl -out ${blastDBName}.blast_db
    """
else
    """
    touch done
    """

}


//Channel.watchPath(kmercountPath+kmercountName, 'create,modify')
genomeChan5.map{it -> file(kmercountPath+kmercountName)}
                    .splitFasta(by: 50000,file: true)
                    .into {kmerFasta; kFast2}
                    //.set { kmerFasta }

kFast2.subscribe {println it}

blast_result = Channel.create()

process BlastOff {

clusterOptions = { writeBlast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 16 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    file 'query.fa' from kmerFasta
    //each genomeName from genomeChan55.toList()

output:
    file blast_result
    //val genomeName into genomeChan6
    //${kmercountPath}${kmercountName}

script:
if(writeBlast == 1)
    """
    #!/bin/bash
    module load blast+/2.6.0
    #cd ${workingDir}
    blastn -db ${workingDir}/${blastDBName}.blast_db -query query.fa -task "blastn-short" -outfmt 6 -num_threads 15 -evalue 1e-2 > blast_result
    """
else
    """
    touch blast_result
    """

}

if (writeBlast == 1){
blast_result.collectFile(name: blastPath + blastName)
            .map {file -> genomeSplitName}
            .set {genomeChan6}
            }
else {
blast_result.collectFile()
            .map {file -> genomeSplitName}
            .set {genomeChan6}
}

genomeChan7 = Channel.create()

process blast2bed {
executor = 'local'
clusterOptions = { b2b == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan6

output:
    val genomeName into genomeChan7


script:
if(b2b == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py blast2bed ${blastPath}${blastName}
    """
else
    """
    touch done
    """

}

genomeChan8 = Channel.create()
//genomeChan85 = Channel.create()

process genClusterMatrix_kmerPrevalence {

clusterOptions = { genMat == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan7

output:
    val genomeName into genomeChan8
    //val genomeName into genomeChan85


script:
if(genMat == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py genClusterKmer ${kmercountPath} ${save} ${genomeName}
    """
else
    """
    touch done
    """

}

reduction_techniques = ['factor','kpca']
//,'feature']

process transform_main {

executor = 'local'

clusterOptions = { trans == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan8
    each technique from reduction_techniques

script:
if(trans == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py transform_main 1 ${reclusterPath} ${technique}
    """
else
    """
    #!/bin/bash
    cd ${workingDir}
    touch *transformed3D.npy
    """

}


//peaks4 = peaks3.map { it -> findPeak(it) }
bestKmerMatrices = Channel.watchPath(reclusterPath+'*.npz','create,modify')
                          .map { file -> findPeak(file.name) }

process transform {

clusterOptions = { trans2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val kmerMat from bestKmerMatrices
    each technique from reduction_techniques

//output:
//    val peakName into peaks3


script:
if(trans2 == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py transform_plot ${kmerMat} ${reclusterPath} ${technique}
    """
else
    """
    #!/bin/bash
    cd ${workingDir}
    touch *transformed3D.npy
    """

}

transformedData = Channel.watchPath('*transformed3D.npy','create,modify')
                         .unique()

process cluster {

clusterOptions = { clust == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    file transformedData

//output:
//    val peakName into peaks3


script:
if(clust == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py cluster ${transformedData} ${reclusterPath} ${best50kmerPath}
    """
else
    """
    touch done
    """

}

