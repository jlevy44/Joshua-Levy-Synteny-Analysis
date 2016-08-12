fileTest=open('testSaveProgram.txt','r')

chunkList=[]
newChunk=[]
grabLine=['$','$']
for line in fileTest:
    if grabLine[0] in line:
        newChunk=[line]
    if grabLine[1] in line:
        newChunk+=[line]
        chunkList+=[newChunk]
        newChunk=[]
        grabLine=['$','$']
    if '[' in line:
        lineList=line.split()
        grabLine = [lineList[2],lineList[4]]

fileTest.seek(0)

findGenes=[]
#Bd2,510683,Bradi2g41430.1,Chr1,2751500,LOC_Os01g38960.1
for chunk in chunkList:
    line1List=chunk[0].split()
    line2List=chunk[1].split()
    line1ImportantInfo=line1List[0].split(',')
    line2ImportantInfo = line2List[0].split(',')
    findGenes+=[[(line1ImportantInfo[0],line1ImportantInfo[2],line2ImportantInfo[2]),
                 (line1ImportantInfo[3],line1ImportantInfo[5],line2ImportantInfo[5])]]

# [[('Bd2', 'Bradi2g41430.1', 'Bradi2g52197.1'), ('Chr1', 'LOC_Os01g38960.1', 'LOC_Os01g57940.1')], [('Bd2', 'Bradi2g52210.1', 'Bradi2g60620.1'), ('Chr1', 'LOC_Os01g58070.1', 'LOC_Os01g71790.1')]]

fileTest.close()

print findGenes