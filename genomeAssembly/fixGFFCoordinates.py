import subprocess, sys
"""python fixgff geneNaming"""
geneNaming = sys.argv[1]
runCommand = lambda x: subprocess.call(x,shell=True)
replaceGFF = 'mv %s %s'%(geneNaming + '.gff',geneNaming + '.gff3')
runCommand('gffread -E %s.gff3 -o- > %s.gff'%(geneNaming,geneNaming))
runCommand(replaceGFF)
with open(geneNaming + '.gff3','r') as f, open(geneNaming + '.gff','w') as f2:
     for line in f:
             if line and len(line.split('\t')) > 3:
                     if '0' == line.split('\t')[3]:
                             f2.write('\t'.join(line.split('\t')[0:3] + ['1'] + line.split('\t')[4:]))
                     else:
                             f2.write(line)

runCommand(replaceGFF)