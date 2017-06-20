def bed2link(bedfile,BedInputPath='./',LinkOutputPath='./'):
    """bed2link inputs list of bed files, their input path, and output path for link file.
    Converts bed files to link.txt files to be used for circos analysis. Files are outputted to specified link out path
    DO NOT USE THIS FOR WHEAT"""
    listLinks = []
    # if hybridum analysis
    listHybridumSubgenomes = ['000']#['001','002','003','004','005','006','007','008','009','010','011','012']

    outputLinkFilename = bedfile.replace('.bed','.link.txt') # link file type
    listLinks.append(outputLinkFilename) # add new file name to list of link files
    open(LinkOutputPath+outputLinkFilename, 'w').close() # create output link file
    # edit link output file
    outputLinkFile = open(LinkOutputPath+outputLinkFilename,'w')
    # input bed file
    inputFile = open(BedInputPath+bedfile,'r')
    # convert each line in each bed file into proper output for link file
    for line in inputFile:
        """ CHANGING:
        523-Chr09N	105496871	106056711	314-Bd1-64279956-64717012
        INTO:
        Chr09N 105496871	106056711 Bd1 64279956 64717012
        """
        lineList = line.split()
        lineList2 = lineList[-1].split('-')
        #if len(lineList2)==5:
        #    outputTuple = (lineList[0].split('_')[1],lineList[1],lineList[2],lineList2[1]+'_'+lineList2[2],lineList2[3],
        #               lineList2[4])
        #else:
        # proper output, see above^^^^

        # for hybridum analysis
        outTupleChr1 = lineList[0].replace('~','-').replace('-','')#[lineList[0].find('-')+1:].replace('~','-').replace('-','')
        outTupleChr2 = lineList2[0]+lineList2[1].replace('~','-')
        for species in listHybridumSubgenomes:
            if species in bedfile.split('-')[0]:
                outTupleChr1 = outTupleChr1.replace('Bh','Bh%s'%species[1:])
            if species in bedfile.split('-')[1]:
                outTupleChr2 = outTupleChr2.replace('Bh', 'Bh%s' % species[1:])


        outputTuple = (outTupleChr1,lineList[1],lineList[2],
                       outTupleChr2,lineList2[2],lineList2[3])
        outputLinkFile.write('%s %s %s %s %s %s\n'%outputTuple)

    inputFile.close()
    outputLinkFile.close()



    """
    PAC4GC.523-PAC2_0.323_5.bed
    PAC4GC.524-PAC2_0.323_5.bed
    """

    """523-Chr09N	105496871	106056711	314-Bd1-64279956-64717012
    523_Chr09N	121481338	121955050	314_Bd1_67267681_67690088
    523_Chr09N	107793319	108299902	314_Bd1_70723158_71434128
    523_Chr07N	61058128	61791054	314_Bd5_18578168_19092045
    523_Chr02N	91071082	91917422	314_Bd1_19591113_20483272 """