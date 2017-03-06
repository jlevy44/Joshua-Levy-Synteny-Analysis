#!/usr/bin/env python
import os 
import sys
import re
from argparse import ArgumentParser
import glob

parser = ArgumentParser()
parser.add_argument("indir")
parser.add_argument("newick")
parser.add_argument("-t", "--time_string", dest = "time_string", default = "96:00:00")
parser.add_argument("-m","--max_threads", dest = "max_threads", type=int, default = 16)
parser.add_argument("-o","--overwrite", dest = "overwrite", action = "store_true", default = False)
options = parser.parse_args()
haldir = os.path.join(options.indir, "hal")
if not os.path.exists(haldir):
   os.makedirs(haldir)

batch_header = """
#$ -l ram.c=5G
#$ -pe pe_slots %s 
#$ -l h_rt=%s
export _JAVA_OPTIONS="-Xmx155g"
cd /global/homes/p/phillips/software/progressiveCactus/
""" % (options.max_threads, options.time_string)

def existCheck(infile):
   return (os.path.exists(infile) and (os.path.getsize(infile) > 32))

seqfiles = glob.glob(os.path.join(options.indir, "*", "seqfile"))
cmd = "bin/runProgressiveCactus.sh --maxThreads %s %s %s %s >& %s" 
if options.overwrite:
   cmd = "bin/runProgressiveCactus.sh --maxThreads --overwrite %s %s %s %s >& %s"
for sf in seqfiles:
   work_dir = os.path.dirname(sf)
   batchdir = os.path.join(work_dir, "batch")
   fname = os.path.basename(work_dir)
   fname = re.sub("FastaOut", "fo", fname)
   if not os.path.exists(batchdir):
      os.makedirs(batchdir)
   batch_file = sf + ".batch.sh"
   err_file = batch_file + ".err"
   out_file = batch_file + ".out"
   stdout_file = batch_file + ".stdout"
   batch_fh = open(batch_file, 'w')
   hal_file = os.path.join(work_dir, os.path.basename(os.path.dirname(sf)) + ".hal")
   if not existCheck(hal_file):
      batch_fh.write(batch_header)
      print >> batch_fh, cmd % (str(int(options.max_threads)), sf, work_dir, hal_file, stdout_file)
      print >> batch_fh, "cp %s %s" % (hal_file, os.path.abspath(haldir))
      batch_fh.close()
      qsub_cmd = "qsub -P plant-analysis.p -N cac_%s -e %s -o %s %s -m abes -M jlevy@lbl.gov" % (fname, err_file, out_file, batch_file) 
      os.system(qsub_cmd)


