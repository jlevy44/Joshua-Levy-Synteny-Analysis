def parseConfigFindPath(stringFind,configFile):
    """findPath will find path or value of associated specified string or info from config file"""
    for line in configFile:
        if stringFind in line: # if find string specified, return pathname or specific value trying to find
            configFile.seek(0)
            return line.split()[-1].strip('\n')
    configFile.seek(0)

with open('kmerCountCompare.txt','r') as f:
    try:
        slurm, local = tuple([int(parseConfigFindPath(x,f)) for x in ['slurm','local']])
    except:
        slurm, local = (0,0)

if slurm == 1:
    executor = "'slurm'"
    clusterOptions = ''#'-N 3 -p regular -D . '
elif (slurm, local) == (0,1):
    executor = 'local'
    clusterOptions = '-P plant-analysis.p -cwd -l exclusive.c -pe pe_slots 3 -e OutputFile.txt'
else:
    executor = 'sge'
    clusterOptions = '-P plant-analysis.p -cwd -l exclusive.c -pe pe_slots 3 -e OutputFile.txt'

with open('nextflow.config','r') as f:
    configLines = f.readlines()

for i, line in enumerate(configLines):
    if 'executor =' in line:
        configLines[i] = 'executor = %s'%(executor)

with open('nextflow.config','w') as f:
    f.write('\n'.join(configLines))