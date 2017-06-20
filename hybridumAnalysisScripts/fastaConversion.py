from os import listdir
from pyfaidx import Fasta, wrap_sequence



fileList = listdir('.')
fileList.remove('BhybridumS_002_v1.0.softmasked.fa')
fileList.remove('q.PAC4GC.hybridum.S.chr.list')
fileList = [fasta for fasta in fileList if fasta.endswith('.fasta')]
with open('q.PAC4GC.hybridum.S.chr.list','r') as f:
    chrList = f.read().split('\n')

print chrList

chrList.remove('')

for file in fileList:
    newfilename = file[:file.find('.')]+'S.fasta'
    print newfilename
    fasta = Fasta(file)
    open(newfilename,'w').close()
    with open(newfilename,'w') as fOut:
        print chrList
        print file,newfilename
        for chromName in chrList:
            print '>'+chromName
            fOut.write('>' + chromName +'\n%s\n'%(''.join(line for line in wrap_sequence(60, str(fasta[chromName])))))