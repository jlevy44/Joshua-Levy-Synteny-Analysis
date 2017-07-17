from pyfaidx import Fasta
from pybedtools import BedTool
import subprocess

badChrM1 = ['3_76388985_88434194']

badChrM2 = ['3_99309338_112736828','3_93379329_144019348']

badInterval = BedTool('\t'.join(['3','76388985','144019348']),from_string=True)

with open('Maize1_SplitInfo.txt','r') as f:
    (BedTool('\n'.join([line[line.find('\t')+1:].strip('\n') + '\t' + line.split('\t')[0] for line in f.readlines()[1:] if line]),from_string=True)-badInterval).sort().saveas('M1.bed')

with open('Maize2_SplitInfo.txt','r') as f:
    (BedTool('\n'.join([line[line.find('\t')+1:].strip('\n') + '\t' + line.split('\t')[0] for line in f.readlines()[1:] if line]),from_string=True)-badInterval).sort().saveas('M2.bed')

subprocess.call('bedtools getfasta -fi Zmays_xxx_Smv4.fa -fo ZmaysM1_640_AGPv4NoAmbiguous.fa -bed M1.bed -name',shell=True)
subprocess.call('bedtools getfasta -fi Zmays_xxx_Smv4.fa -fo ZmaysM2_641_AGPv4NoAmbiguous.fa -bed M2.bed -name',shell=True)