def bed2link(listBedFiles,BedInputPath,LinkOutputPath):
    listLinks = []
    for file in listBedFiles:
        outputLinkFilename = file.replace('.bed','.link.txt')
        listLinks.append(outputLinkFilename)
        open(LinkOutputPath+outputLinkFilename, 'w').close()
        outputLinkFile = open(LinkOutputPath+outputLinkFilename,'w')
        inputFile = open(BedInputPath+file,'r')
        for line in inputFile:
            lineList = line.split()
            lineList2 = lineList[-1].split('-')
            #if len(lineList2)==5:
            #    outputTuple = (lineList[0].split('_')[1],lineList[1],lineList[2],lineList2[1]+'_'+lineList2[2],lineList2[3],
            #               lineList2[4])
            #else:
            outputTuple = (lineList[0].split('-')[1],lineList[1],lineList[2],lineList2[1],lineList2[2],lineList2[3])
            outputLinkFile.write('%s %s %s %s %s %s\n'%outputTuple)

        inputFile.close()
        outputLinkFile.close()

    return listLinks


    """
    PAC4GC.523-PAC2_0.323_5.bed
    PAC4GC.524-PAC2_0.323_5.bed
    """

    """523-Chr09N	105496871	106056711	314-Bd1-64279956-64717012
    523_Chr09N	121481338	121955050	314_Bd1_67267681_67690088
    523_Chr09N	107793319	108299902	314_Bd1_70723158_71434128
    523_Chr07N	61058128	61791054	314_Bd5_18578168_19092045
    523_Chr02N	91071082	91917422	314_Bd1_19591113_20483272 """