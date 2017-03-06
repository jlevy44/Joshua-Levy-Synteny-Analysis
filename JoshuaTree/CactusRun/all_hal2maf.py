#!/usr/bin/env python
import os
import sys
import glob
import multiprocessing as mp
from argparse import ArgumentParser
from time import sleep

def main():
   parser = ArgumentParser()
   parser.add_argument("indir")
   parser.add_argument("outdir")
   parser.add_argument("hal2maf_path", default="~/software/progressiveCactus/submodules/hal/bin/hal2maf")
   parser.add_argument("-n", "--numcpu", dest="numcpu", default=12, type=int)
   options = parser.parse_args()

   if not os.path.exists(options.outdir):
      os.makedirs(options.outdir)

   for halfile in glob.glob(os.path.join(options.indir, "*.hal")):
      print halfile
      proc = mp.Process(target=runhal, args = ([halfile, options]))   
      proc.daemon = True
      proc.start()
      while len(mp.active_children()) > options.numcpu:
         sleep(1)
   while len(mp.active_children()) > 0:
      sleep(1)

def runhal(halfile, options):
   bname = os.path.splitext(os.path.basename(halfile))[0]
   outfile = os.path.join(options.outdir, bname + ".maf")
   os.system("%s %s %s" % (options.hal2maf_path, halfile, outfile))

if __name__ == '__main__':
       main()
