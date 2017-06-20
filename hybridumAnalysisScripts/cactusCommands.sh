#!/bin/bash
python syn_fasta_to_cactus.py /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/fasta/ /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/ grassnewick.nh -p prot_dict
mkdir /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/hal
python runcac.py -t 11:59:30 /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/ /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/root/grassnewick.nh
wait
mkdir /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/maf
python all_hal2maf.py /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/hal /Users/jlevy/Desktop/Projects/Synteny/OtherProjects/output/maf -t 12