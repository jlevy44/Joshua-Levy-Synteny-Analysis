"""Files are in: /projectb/sandbox/plant/synfasta_scripts

Step 1: Convert fastas from Josh's output to cactus-ready directories. Example:
python syn_fasta_to_cactus.py /global/dna/projectdirs/plant/pangenomics/CNS/sean_hybridum_final_synteny/FastaFileOutputsSsubgenome/ /projectb/sandbox/plant/synfasta_scripts/output grassnewick.nh -p prot_dict
(better to use a more informative directory name than my choice of "output")

The prot_dict file contains the number to name mapping for any proteome that is not in PAC, edit as needed.

Step 2: Queue up cactus jobs:
python runcac.py -t 11:59:30 /projectb/sandbox/plant/synfasta_scripts/output/ /projectb/sandbox/plant/synfasta_scripts/grassnewick.nh

Try running with that time string (11:59:30) to stay on the short queue first and see what will finish in 12 hours. Then run again with 72 or 96 hours (go higher if you expect very long synteny blocks). The highest you can go is -t 240:00:00 for the max of 10 days. If you need to do that, expect your jobs to be in the queue for a very long time before they start running.

The finished hal files will be copied over to the "hal" directory wherever your cactus dirs are, so  /projectb/sandbox/plant/synfasta_scripts/output/hal, in this case.

Step 3: Convert to maf if you want
python all_hal2maf.py <directory_with_halfiles> <maf_output_directory> -t 12
The -t is optional; this will run on a gpint multithreaded, default is to use 12 cores, which will finish in an hour or so."""
import subprocess, os, shutil, sys

def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of relevant filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of files or list information
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    return listOfItems

def parseConfigFindPath(stringFind,configFile):
    """findPath will find path or value of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or specific value trying to find
            configFile.seek(0)
            return line.split()[-1].strip('\n')

config = open('/Users/jlevy/Desktop/Projects/Synteny/OtherProjects/cactusConfig.txt','r')
root = parseConfigFindPath('root',config)
fastaFolder = parseConfigFindPath('fastaPath',config)
outputPath = parseConfigFindPath('outputPath',config)
config.close()
open('cactusCommands.sh','w').close()
Commands = ['#!/bin/bash\n','python syn_fasta_to_cactus.py %s %s grassnewick.nh -p prot_dict\n'%(fastaFolder,outputPath),
            'wait\n','mkdir %shal\n'%outputPath,'python runcac.py -t 11:59:30 %s %sgrassnewick.nh\n'%(outputPath,root),'wait\n',
            'mkdir %smaf\n'%outputPath,'python all_hal2maf.py %shal %smaf -t 12'%(outputPath,outputPath)]
cactusComm = open('cactusCommands.sh','w')
cactusComm.writelines(Commands)
cactusComm.close()
subprocess.call('nohup sh cactusCommands.sh', shell=True)
