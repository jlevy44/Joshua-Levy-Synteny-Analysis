#!/bin/bash
module unload gcc 
module load gcc/5.4.0
nohup python /projectb/sandbox/plant/synfasta_scripts/all_hal2maf.py output/hal output/maf /global/homes/p/phillips/software/progressiveCactus/submodules/hal/bin/hal2maf -n 12
