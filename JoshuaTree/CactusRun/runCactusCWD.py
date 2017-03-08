import os,subprocess

"""
python /global/dna/projectdirs/plant/pangenomics/synfasta_scripts/jlphillips_lbl-synfasta_scripts-31de175504f3/syn_fasta_to_cactus.py -n 5 -p prot_dict /global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/CactusRun/FastaFiles /global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/CactusRun/output /global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/CactusRun/grassnewick.nh


python runcac.py -t 138:59:30 /global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/CactusRun/output/ /global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/CactusRun/grassnewick.nh
"""

syn2fasta = '/global/dna/projectdirs/plant/pangenomics/synfasta_scripts/jlphillips_lbl-synfasta_scripts-31de175504f3/syn_fasta_to_cactus.py'
runcac = '/global/dna/projectdirs/plant/pangenomics/synfasta_scripts/jlphillips_lbl-synfasta_scripts-31de175504f3/runcac.py'
cwd = os.getcwd()+'/'

subprocess.call([syn2fasta,'-n','5','-p','prot_dict',cwd+'FastaFiles',cwd+'output',cwd+'grassnewick.nh'])
subprocess.call([runcac,'-t','138:59:30',cwd+'output',cwd+'grassnewick.nh'])
