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
constructSample = findValue('constructSample ').asType(Integer);

chanBuildSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
//.toList()
//chanBuildSamples.subscribe { println "value: $it" }
workingDir = new File('').getAbsolutePath()
chanBuildRef = Channel.fromPath('referenceGenomes' + '/*',type: 'dir', relative: true)
//= Channel
                    //
                    /* .subscribe(onNext: { println it}, onComplete: { println 'Done.' }) */



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
python ${workingDir}/writeShFiles.py
"""
else {
    """
    #!/bin/bash
    echo Not writing analysis files...
    touch done
    """
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

clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6 -e OutputFile.txt'

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




linkageChannel1 = linkageChannel1.take(1)


linkageChannel2 = Channel.create()

process buildSample{
clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6 -e OutputFile.txt'

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
"""touch done"""
}
//linkageChannel2 = linkageChannel2.take(1)


process nucmerfy {

clusterOptions = '-P plant-analysis.p -cwd -l high.c -pe pe_slots 6 -e OutputFile.txt'

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
if(constructSample)
"""
#!/bin/bash
touch done
cd ${workingDir}/${version}/${sample}
sh nucCommand.sh
"""
else
"""touch done"""
}

linkageChannel3 = Channel.create()
linkageChannel4 = Channel.create()
linkageChannel5 = Channel.create()
refStrList = refStr.toList()
//linkageChannel3 = linkageChannel3.take(1)
process com_1_2 {

input:
val sample from linkageChannel3
each refText from refStrList

//file 'done' from linkageChannel3

//each sample from chanBuildSamples

output:
val sample into linkageChannel4
//file 'done' into linkageChannel4


script:
if(constructSample)
"""
#!/usr/bin/python
import os, subprocess, sys
open('done','w').close()
os.chdir('../../..')
from pipelineFunctions import *
os.chdir(${version}/${sample})
if all([os.path.isfile('%s.%s.lifted.anchors' %(${sample}, ref)) and os.stat('%s.%s.lifted.anchors' %(${sample}, ref)).st_size > 0 for ref in ${refText}]) == 0:
    subprocess.call('sh constructv1_1.sh',shell=True)
    sampleCount = 0
    for ref in weights.keys():
        sampleCount = replaceGeneNames(${sample}, ref, sampleCount)
try:
    tiling2bed('%snuc.tiling'%${CDS}, ${CDS}, ${sample}, ${sample}+'_%ssyn'%${CDS}+'.bed')
except:
    print sys.exc_info()[0]
replaceGeneNames(${sample},${CDS},0,1)
subprocess.call('sh constructv1_2.sh',shell=True)
"""
else
"""touch done"""
}
//linkageChannel4 = linkageChannel4.take(1)
process allmaps {
clusterOptions = '-P plant-analysis.p -cwd -l h_rt=50:00:00 -pe pe_slots 32 -e OutputFile.txt'
queue 'long'

input:
val sample from linkageChannel4
//each sample from chanBuildSamples
//file 'done' from linkageChannel4


output:
val sample into linkageChannel5



script:
//if(File('${version}/${sample}/multipleMapping.bed').length() == 0)
if(constructSample)
"""
#!/bin/bash
touch done
cd ../../../${version}/${sample}
sh qsub_build.sh
"""
else
    error "Unable to build new genome: ${sample}"
}
