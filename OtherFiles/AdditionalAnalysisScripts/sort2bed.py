# convert sort2 files to bed files
listSort2Files = ['q.PAC4GC.523.sort2','q.PAC4GC.524.sort2']

for file in listSort2Files:
    sortReadFile = open(file,'r')
    open(file.replace('sort2','bed'),'w').close()
    bedOutFile = open(file.replace('sort2','bed'),'w')
    for line in sortReadFile:
        lineList = line.split()
        lineList[1] = file.split('.')[-2]+'_'+lineList[1]
        bedOutFile.write('%s\t%s\t%s\n'%(tuple(lineList[1:])))

    sortReadFile.close()
    bedOutFile.close()