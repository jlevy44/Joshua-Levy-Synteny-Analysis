import sys

inputList = sys.argv

for e in inputList:
    if e.endswith('.gff3') or e.endswith('.gff'):
        outfilename = e.replace('.gff3','.CDS.bed').replace('.gff','.CDS.bed')
        with open(e,'r') as f:
            #for line in f:
            #    if line and '#' not in line and 'CDS' == line.split()[2]:
            #        print line.split('\t')[8].split(';')[1].replace('Name=','')
                    #print (line.split()[0],line.split()[3],line.split()[4],line.split()[8].split(';')[1].replace('Name=',''))
            #4719-83673	Bdist012v0	CDS	6963	7181	100	+	0	ID=Bdist012v02g37802.1.mrna1.cds1;Name=Bdist012v02g37802.1;Parent=Bdist012v02g37802.1.mrna1;Target=Bdist012v02g37802.1 1 219 +
            outTuples = [(line.split('\t')[0],line.split('\t')[3],line.split('\t')[4],line.split('\t')[8].split(';')[1].replace('Name=','')) for line in f if line and '#' not in line and 'CDS' == line.split('\t')[2]]#.split(';')[0].replace('Parent=','').strip('\n'),'gene='+line.split('\t')[8].split(';')[0].replace('Parent=','').strip('\n').replace('.mrna1','').replace('.mrna2','')) for line in f if line and '#' not in line and 'CDS' == line.split('\t')[2]]
        open(outfilename,'w').close()
        #print outTuples
        with open(outfilename,'w') as f:
            f.write('\n'.join('\t'.join(part for part in outtuple) for outtuple in outTuples))