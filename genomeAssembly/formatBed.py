import sys, os
from pybedtools import BedTool

refSamp = sys.argv[1] # r or s

fname = sys.argv[2]+'.bed'

version = sys.argv[3] # v0, v1, ...

if refSamp == 'r':
    folder = './referenceGenomes/'
else:
    folder = './'+version+'/'+sys.argv[2]+'/'
try:
    sys.argv[4] == '1'
except:
    os.chdir(folder)

with open(fname,'r') as f:
    bedLines = f.readlines()

with open(fname,'w') as f:
    for line in bedLines:
        if line and ',' not in line.split('\t')[3]:
            if '.' in line:
                f.write(line[:line.find('\t')]+line[line.find('\t'):line.find('.')].replace('-','')+'\t'+line[line.rfind('0'):].strip('\n') + '\n')
            else:
                f.write(line.replace(line.split('\t')[3],line.split('\t')[3].replace('-',''))+'\n')

if refSamp == 'r':
    BedTool(fname).sort().merge(c=[4,5,6],o=['distinct','distinct','distinct']).saveas(fname)
elif refSamp == 's':
    b = BedTool(fname).sort().merge(c=4,o='count_distinct')
    a=''
    for line in str(b).split('\n'):
        if line and int(line.split('\t')[-1]) > 1:
            a += line + '\n'
    BedTool(fname).subtract(BedTool(a,from_string=True)).sort().merge(c=[4,5,6],o=['distinct','distinct','distinct']).saveas(fname)

with open(fname,'r') as f:
    bedLines = f.readlines()

with open(fname,'w') as f:
    for line in bedLines:
        if line and ',' not in line.split('\t')[3]:
            f.write(line)

