#!/usr/bin/env nextflow

String findValue(String inputStr) {
configFile = file('scaffoldConfig.txt')
allLines = configFile.readLines()
for( line in allLines ) {
    if(line.startsWith(inputStr)){
        println (line - inputStr - '= ')
        return (line - inputStr - '= ');
    }
}

}

params.folder = './'
//String workingDir = System.getProperty("user.dir");
writeSh = findValue('writeSh ').asType(Integer);
println 'writeSh ' + writeSh.asType(String)
buildRef = findValue('buildRef ').asType(Integer);
version = findValue('version ' );
CDS = findValue('CDS ' );
CDSFasta = findValue('CDSFasta ' );
geneNameOld = findValue('geneNameOld ');
buildSamp = findValue('buildSample ').asType(Integer);
BB = findValue('BB ').asType(Integer);
nuc = findValue('nuc ').asType(Integer);
com1_2 = findValue('com1_2 ').asType(Integer);
allmaps = findValue('allmaps ').asType(Integer);
nextVersion = 'v' + ((version - 'v').asType(Integer)+1).asType(String)
chanBuildSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
//.toList()
//chanBuildSamples.subscribe { println "value: $it" }
workingDir = new File('').getAbsolutePath()
nextDir = file(workingDir+'/'+nextVersion)
result = nextDir.mkdir()
chanBuildRef = Channel.fromPath('referenceGenomes' + '/*',type: 'dir', relative: true)
chanNextSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
                        .map { nextVersion + '/' + it.name - version + nextVersion}
                        .map{ file(it) }
                        .subscribe { it.mkdirs() }
chanNextSamples2 = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
                        .map { nextVersion + '/' + it.name - version + nextVersion}
                        .map{ file(it) }
                        .subscribe { println it }
//= Channel
                    //
                    /* .subscribe(onNext: { println it}, onComplete: { println 'Done.' }) */

//bridgeChannel = Channel.create()
/*
process newVersionFolders {
executor = 'local'

input:
    file sample from chanNextSamples

output:
    file 'done' into bridgeChannel

script:
"""
#!/bin/bash
touch done
mkdir ${workingDir}/${nextVersion}/${sample}
"""
}
*/

linkageChannel = Channel.create()
process writeShFiles {

executor = 'local'

output:
    file 'done' into linkageChannel

script:
if(writeSh)
"""
#!/bin/bash
touch done
module load bedtools/2.25.0
python ${workingDir}/writeShFiles.py
cd ${workingDir}/referenceGenomes/${CDS}
python ../../generateCentromereBeds.py
cd -
"""
else {
    """
    #!/bin/bash
    echo Not writing analysis files...
    touch done
    """
    //mkdir ${workingDir}/${version}
    }
}


linkageChannel = linkageChannel.take(1)
linkageChannel15= Channel.create()
refStr = Channel.create()

process findRef {
executor = 'local'
input:
    file 'done' from linkageChannel

output:
    file 'done' into linkageChannel15
    stdout refStr

script:
"""
#!/usr/bin/python
open('touch','w').close()
with open('../../../references.txt','r') as f:
    print f.read()
"""
/*
refFile = file('../../../references.txt')
allLines = refFile.readLines()
for( line in allLines ) {
    if(line) {
    refText = line
    break;
    }
}

println refText
done = file('\${NXF_WORK}/done')
*/
}


//refText = "[${chanBuildRef.join(',')}]"
//println refText


linkageChannel15 = linkageChannel15.take(1)



linkageChannel1 = Channel.create()

process buildReference {

clusterOptions = { buildRef == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
file 'done' from linkageChannel15
each ref from chanBuildRef.toList()

output:
file 'done' into linkageChannel1


script:
if(buildRef)
"""
#!/bin/bash
touch done
cd ${workingDir}/referenceGenomes/${ref}
pwd
sh buildRef.sh
"""
else
"""touch done"""
}




linkageChannel1 = linkageChannel1.last()


linkageChannel2 = Channel.create()

process buildSample{
clusterOptions = {buildSamp == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
each sample from chanBuildSamples.toList()
file 'done' from linkageChannel1

output:
val sample into linkageChannel2
//file 'done' into linkageChannel2

//when:
//Channel.fromPath('${version}/${sample}/*${sample}.cds').ifEmpty('E').toList() == ['E'] || Channel.fromPath('${version}/${sample}/*${sample}.bed').ifEmpty('E').toList()
//file('${version}/${sample}/${sample}.cds').exists() == 0 || file('${version}/${sample}/${sample}.bed').exists() == 0

script:
if(buildSamp)
"""
#!/bin/bash
touch done
cd ${workingDir}/${version}/${sample}
pwd
sh build.sh
"""
else
"""#!/bin/bash
touch done"""
}
//linkageChannel2 = linkageChannel2.take(1)

linkageChannel3 = Channel.create()
process nucmerfy {

clusterOptions = {nuc == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}
//-l exclusive.c
input:
val sample from linkageChannel2
//file 'done' from linkageChannel2
//each sample from chanBuildSamples
//val chanBuildSamples
//file 'done'


output:
val sample into linkageChannel3
//file 'done' into linkageChannel3
//file 'done'

//when:
//Channel.fromPath('${version}/${sample}/*${CDS}nuc.tiling').ifEmpty('E').toList() == ['E']
//file('${version}/${sample}/*${CDS}nuc.tiling').exists() == 0

script:
if(nuc)
"""
#!/bin/bash
touch done
cd ${workingDir}/${version}/${sample}
sh nucCommand.sh
"""
else
"""touch done"""
}


linkageChannel3_5 = Channel.create()
process BBfy {

clusterOptions = {BB == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 8 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}
//-pe pe_slots 6
input:
val sample from linkageChannel3



output:
val sample into linkageChannel3_5


script:
if(BB)
"""
#!/bin/bash
export _JAVA_OPTIONS='-Xmx15G'
touch done
module load bedtools/2.25.0
cd ${workingDir}/${version}/${sample}
sh BB_build.sh > nohup.out
"""
else
"""touch done"""
}


linkageChannel4 = Channel.create()
linkageChannel5 = Channel.create()
refStrList = refStr.take(1)
    .map { it - '\n' }
    .toList()
//linkageChannel3 = linkageChannel3.take(1)
process com_1_2 {
clusterOptions = { com1_2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 3 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}
//module = 'bedtools/2.25.0'

input:
val sample from linkageChannel3_5
each refText from refStrList

//file 'done' from linkageChannel3

//each sample from chanBuildSamples

output:
val sample into linkageChannel4
//file 'done' into linkageChannel4


script:
if(com1_2)
"""
#!/bin/bash
touch done
module load bedtools/2.25.0
cd ${workingDir}
python com1_2.py ${sample} ${refText} ${CDS} ${version}
"""
else
"""touch done"""
}
//linkageChannel4 = linkageChannel4.take(1)
process Allmaps {
clusterOptions = '-P plant-analysis.p -cwd -l h_rt=24:00:00 -pe pe_slots 16 -e OutputFile.txt'
queue = 'long'
errorStrategy = 'ignore'

input:
val sample from linkageChannel4
//each sample from chanBuildSamples
//file 'done' from linkageChannel4


output:
val sample into linkageChannel5



script:
//if(File('${version}/${sample}/multipleMapping.bed').length() == 0)
if(allmaps)
"""
#!/bin/bash
touch done
cd ../../../${version}/${sample}
sh qsub_build.sh
"""
else
    error "Unable to build new genome: ${sample}"
}
