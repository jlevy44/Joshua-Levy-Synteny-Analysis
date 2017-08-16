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
clust = findValue('ClusterAll ').asType(Integer);


genomeSplitName = genome - '_split' - '.fa' + '_split.fa';
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
    val genomeName into genomeChan55


script:
if(fromFasta == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    python subgenomeClusteringInterface.py kmer2Fasta ${kmercountPath}
    makeblastdb -in ${genomeFullPath} -dbtype nucl -out ${genomeName}.blast_db
    """
else
    """
    touch done
    """

}


//Channel.watchPath(kmercountPath+kmercountName, 'create,modify')
genomeChan5.map{it -> file(kmercountPath+kmercountName)}
                    .splitFasta(by: 50000,file: true)
                    .set { kmerFasta }

blast_result = Channel.create()

process BlastOff {

clusterOptions = { writeBlast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    file 'query.fa' from kmerFasta
    val genomeName from genomeChan55

output:
    file blast_result
    //val genomeName into genomeChan6
    //${kmercountPath}${kmercountName}

script:
if(writeBlast == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    blastn -db ./${genomeName}.blast_db -query query.fa -task "blastn-short" -outfmt 6 -num_threads 8 -evalue 1e-2 > blast_result
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

process genClusterMatrix_kmerPrevalence {

clusterOptions = { genMat == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan7

output:
    val genomeName into genomeChan8


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

peaks = Channel.create()//genomeChan8.last()
        //            .subscribe {println it}
        //            .fromPath()
                    //.watchPath(peakPath+'*.fa')

process kmerHist {

clusterOptions = { kHist == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val genomeName from genomeChan8

output:
    file 'PeaksOutNames.txt' into peaks


script:
if(kHist == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py kmerRelatedHistogram ${kmercountName} ${save}
    touch PeaksOutNames.txt
    """
else
    """
    #!/bin/bash
    cd ${workingDir}
    touch PeaksOutNames.txt
    """

}

//peaksWait = peaks.last()

peaksFinal = peaks.splitText()
                  .flatMap{it -> tuple(it, it - '.fa' + '.BLASTtsv.txt')}
                //.flatMap{file -> tuple(file.name, file.name - '.fa' + '.BLASTtsv.txt')}

                    //.last()
                    //.subscribe {println it}
                    //.fromPath(peakPath+'*.fa')
                    //.watchPath(peakPath+'*.fa','create')
                    //.flatMap{file -> tuple(file, file.name, file.name - '.fa' + '.BLASTtsv.txt')}
                    //.splitText(by: 5000, file: True)
                    //.set { chunksChannel }

peaks2 = Channel.create()

process blastPeaks {

clusterOptions = { bPeak == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    set peakName, peakBlastName from peaksFinal

output:
    set peakName, peakBlastName into peaks2

//-out ${blastPath}/${peakBlastName} >>
script:
if(bPeak == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    blastn -db ./${genomeName}.blast_db -query ${peakName} -task "blastn-short" -outfmt 6 -num_threads 8 -evalue 1e-2 -out ${blastPath}/${peakBlastName}
    """
else
    """
    touch done
    """

}

// can try .collectFile().subscribe { merged_file -> merged_file.copyTo(out_dir) }



peaks3 = Channel.create()

process genPeakCluster {

clusterOptions = { bPeak2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    set peakName, peakBlastName from peaks2

output:
    val peakName into peaks3


script:
if(bPeak2 == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py peakClusteringMatrix ${kmercountPath} ${peakName} ${peakBlastName} ${save}
    """
else
    """
    touch done
    """

}



peaks4 = peaks3.map { it -> findPeak(it) }

process transform {

clusterOptions = { trans == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
    val peak from peaks4

//output:
//    val peakName into peaks3


script:
if(trans == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py transform_plot ${peak}
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
    python subgenomeClusteringInterface.py cluster ${transformedData}
    """
else
    """
    touch done
    """

}

