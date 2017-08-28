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

    //String[] parts = peakfname.split("_");
    //String peakfinalname = parts[0];
    int finalPosition = peakfname.lastIndexOf('_');
    String peakfinalname = peakfname.substring(0,finalPosition);
    return peakfinalname;

}

// parse config

pythonPath = findValue('pythonPath ' );
systemPath = findValue('systemPath ' );
blastPath = findValue('blastPath ' );
kmercountPath = findValue('kmercountPath ' );
reclusterPath = findValue('reclusterPath ');
best500kmerPath = findValue('kmer500BestPath ');
fastaPath = findValue('fastaPath ' );
bedPath = findValue('bedPath ' );
sortPath = findValue('sortPath ' );
sortedbedPath = findValue('sortedbedPath ' );
genome = findValue('genome ' );
peakPath = findValue(' peakPath');


//preprocess = findValue('preprocess ').asType(Integer);


save = findValue('save ');
old = findValue('old ');

slurm = findValue('slurm ').asType(Integer);

splitFast = findValue('splitFasta ').asType(Integer);
writeKmer = findValue('writeKmer ').asType(Integer);
fromFasta = findValue('kmer2Fasta ').asType(Integer);
original = findValue('original ').asType(Integer);
writeBlast = findValue('writeBlast ').asType(Integer);
b2b = findValue('blast2bed ').asType(Integer);
genMat = findValue('generateClusteringMatrix ').asType(Integer);
kHist = findValue('kmerHistogram ').asType(Integer);
bPeak = findValue('blastPeaks ').asType(Integer);
bPeak2 = findValue('peakClusterMatrix ').asType(Integer);
trans = findValue('transformData ').asType(Integer);
trans2 = findValue('transformDataChildren ').asType(Integer);
clust = findValue('ClusterAll ').asType(Integer);
extract = findValue('extract ').asType(Integer);
kmerBlast = findValue('kmerBlast ').asType(Integer);
kmerGraph = findValue('kmerGraph ').asType(Integer);



genomeSplitName = genome - '_split' - '.fa' + '_split.fa';
blastDBName = genomeSplitName - '.fa'
genomeFullPath = fastaPath + genomeSplitName;
originalGenome = { original ? fastaPath + genome : genomeFullPath };
kmercountName = genomeSplitName - '.fa' + '.kcount' + '.fa';
blastName = kmercountName - '.fa' + '.BLASTtsv.txt';
workingDir = new File('').getAbsolutePath();
blastDBName2 = {original == 1 ? genome - '.fa' : blastDBName }

check = Channel.from(genomeSplitName,genomeFullPath,kmercountName,blastName,workingDir)
                .subscribe{println it}

//gCh = Channel.from(genome)
//            .subscribe {println it + 'aaaa'}

genomeChan = Channel.from(genome - '\n')

genomeChan2 = Channel.create()

process splitFastaProcess {
executor = 'local'
cpus = { splitFast == 1 ? 2 : 1 }
//'-N 2 -p regular -D . '}

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

//clusterOptions = { slurm == 0 ? { writeKmer == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 9 -p regular -D . '}
cpus = { writeKmer == 1 ? 9 : 1 }

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
//clusterOptions = { slurm == 0 ? { fromFasta == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
//clusterOptions = { slurm == 1 ? '-D .' : '-P plant-analysis.p -cwd' }
cpus = { fromFasta == 1 ? 2 : 1 }

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
    makeblastdb -in ${genomeFullPath} -dbtype nucl -out ${blastDBName}.blast_db
    """
else
    """
    touch done
    """

}

process createOrigDB {
executor = 'local'
//clusterOptions = { slurm == 0 ? { original == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { original == 1 ? 2 : 1 }

input:
    val genomeName from genomeChan55

script:
if(original == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    module load blast+/2.6.0
    makeblastdb -in ${originalGenome} -dbtype nucl -out ${blastDBName2}.blast_db
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

//clusterOptions = {slurm == 0 ? { writeBlast == 1 ? '-P plant-analysis.p -cwd -l h_rt=24:00:00 -pe pe_slots 16 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : { writeBlast == 1 ? '-N 16 -p regular -D . ' : '-N 1 -p regular -D . ' }}
cpus = { writeBlast == 1 ? 16 : 1 }


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
//clusterOptions = {slurm == 0 ? { b2b == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { b2b == 1 ? 2 : 1 }


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

//clusterOptions = {slurm == 0 ? { genMat == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { genMat == 1 ? 2 : 1 }


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

/*
if(trans == 0){
genomeChan8 = genomeChan8.take(1)
}
*/

process transform_main {

executor = 'local'

//clusterOptions = { slurm == 0 ? { trans == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { trans == 1 ? 2 : 1 }


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
    touch main_${technique}_transformed3D.npy
    """

}


//peaks4 = peaks3.map { it -> findPeak(it) }
bestKmerMatrices = Channel.watchPath(reclusterPath+'*.npz','create,modify')
                          .map { file -> findPeak(file.name) }
                          .unique()

process transform {

//clusterOptions = { slurm == 0 ? { trans2 == 1 ? '-P plant-analysis.p -cwd -l high.c -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}

cpus = { trans2 == 1 ? 2 : 1 }

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
    touch ${kmerMat}_${technique}_transformed3D.npy
    """

}

transformedData = Channel.watchPath('*transformed3D.npy','create,modify')
                         .unique()
                         .map {file -> file.name - 'transformed3D.npy'}
                         //.filter((file.name).startsWith('main'))

clusterModels = ['KMeans']//,'SpectralClustering']
kmerBest500Files = Channel.create()
process cluster {

//clusterOptions = { slurm == 0 ? { clust == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { clust == 1 ? 2 : 1 }


input:
    val transformedData
    each model from clusterModels

output:
    file 'test.txt' into kmerBest500Files//stdout into kmerBest500Files//val 'kmer500Best_${model}${transformedData}n3.fa' into kmerBest500Files
//${best500kmerPath}/

script:
if(clust == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py cluster ${transformedData}transformed3D.npy ${reclusterPath} ${best500kmerPath} ${model}
    cd -
    echo kmer500Best_${model}${transformedData}n3.fa > test.txt
    """
else
    """
    echo kmer500Best_${model}${transformedData}n3.fa > test.txt
    """

}

subgenomeFolders = Channel.watchPath('analysisOutputs/*.txt')
                          .map {file -> file.name - '.txt'}

process subgenomeExtraction {


//clusterOptions = { slurm == 0 ? { extract == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 9 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 9 -p regular -D . '}
cpus = { extract == 1 ? 9 : 1 }


input:
    val subgenomeFolder from subgenomeFolders

//output:
//    val peakName into peaks3


script:
if(extract == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py subgenomeExtraction ./analysisOutputs/${subgenomeFolder} ${fastaPath} ${genomeSplitName} ${genome}
    """
else
    """
    touch done
    """
}

//kmerBest500Files = Channel.watchPath(best500kmerPath + '*.fa','create,modify')
//                          .unique()
//                         .flatMap { file -> tuple(file.name, file.name - '.fa') }

kmerBest500Files.splitText()
                .filter {it.toString().size() > 1}
                .into {best_kmer;printkmer}

          //.filter {( it =~ /recluster/ ) == 0 &&  it.isEmpty() == 0}

best_kmer.map { it -> tuple(it - '\n', it - '\n' - '.fa') }//.filter({ it })//( it =~ /recluster/ ) == 0 &&  it.isEmpty() == 0})//it.contains('recluster') == 0 && it.isEmpty() == 0})
         .set {best_kmer2}
                //.into {best_kmer2; printFlat}
                //

printkmer.subscribe {println it}

//printFlat.subscribe {println it}

kmer_blasted = Channel.create()

process kmerBlastOff {

//clusterOptions = { slurm == 0 ? { kmerBlast == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 16 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 16 -p regular -D . '}
cpus = { kmerBlast == 1 ? 16 : 1 }


input:
    set query, queryFolder from best_kmer2 // kmerBest500Files

output:
    val queryFolder into kmer_blasted

script:
if(kmerBlast == 1)
    """
    #!/bin/bash
    module load blast+/2.6.0
    cd ${workingDir}
    mkdir ${best500kmerPath}/${queryFolder}
    blastn -db ${workingDir}/${blastDBName}.blast_db -query ${best500kmerPath}/${query} -task "blastn-short" -outfmt 6 -num_threads 15 -evalue 1e-2 > ${best500kmerPath}/${queryFolder}/${queryFolder}.blast.txt
    """
else
    """
    touch blast_result
    """

}

process kmerGraphs {

//clusterOptions = { slurm == 0 ? { kmerGraph == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 2 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' } : '-N 2 -p regular -D . '}
cpus = { kmerGraph == 1 ? 2 : 1 }

input:
    val queryFolder from kmer_blasted

script:
if(kmerGraph == 1)
    """
    #!/bin/bash
    cd ${workingDir}
    python subgenomeClusteringInterface.py generateKmerGraph ${best500kmerPath} ${queryFolder}
    """
else
    """
    touch blast_result
    """

}