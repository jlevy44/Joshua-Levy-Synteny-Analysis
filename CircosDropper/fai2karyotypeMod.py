from random import randint
def fai2karyotype(inputfilename,BPsThreshold=100000):
    """fai2karyotype imports list of faiFiles (convertFiles) and changes them to a list of KaryoteFiles and exports
    karyotype .txt filetypes corresponding to each species to the karyotype file path. For use with circos analysis.

    NOTE: ONLY USE ONE QUERY SPECIES ELSE WILL DAMAGE ANALYSIS. Accepts certain threshold of basepairs for chr/scaffd"""

    # if hybridum
    hybridumAnalyze = 0
    speciesNameId = inputfilename[:inputfilename.rfind('_')]
    protId = speciesNameId.split('_')[-1]
    if speciesNameId.split('_')[0].endswith('D'):
        hybridumAnalyze = 1
        replaceChr = ''
    elif speciesNameId.split('_')[0].endswith('S'):
        hybridumAnalyze = 1
        replaceChr = ''
    #import fai file
    inputFile = open(inputfilename,'r')
    #output karyotype .txt file name
    outputFilename = 'karyotype.'+protId+'.txt'
    # panicum virgatum and subgenomes # better if P. vir is query species, FIXME may need to change if target speci
    #create karyotype file for species; fai file will be converted
    with open(outputFilename,'w') as outputFile:
        # count number of chromosomes/scaffolds (helps determine color of chromosomes and scaffolds)
        chrCount = 1
        scaffoldOutList = []
        # for hybridum Analysis
        for line in inputFile: #FIXME change naming conventions!!!
            # generate output tuple from parsing input from fai
            # Bd2	59130575	76056752	80	81
            # TO
            # chr - hs10 10 0 135037215 chr10
            # For example
            lineList = line.split()
            # for example, consider line with Bd2
            if int(lineList[1])>BPsThreshold: # if length chromosome over certain number BPs specified in config
                outTup1 = lineList[0]
                if hybridumAnalyze:
                    outTup1 = outTup1.replace('Bh','Bh%s'%replaceChr)
                if '_' in lineList[0]:
                    # output name that takes into account scaffolds
                    outputTuple = (protId+outTup1,outTup1[0]+outTup1.split('_')[-1],lineList[1],chrCount)
                    with open('karyotype.'+protId+'.correspondence'+'.txt','a') as f:
                        f.write('original='+outTup1+'\tmodifiedName=%s\n'%(outTup1[0]+outTup1.split('_')[-1]))
                else:
                    # output chromosome name, chr name, how large chromosome is, and color of chromosome indicated by %d
                    outputTuple = (protId+outTup1, outTup1, lineList[1], chrCount)
                if 'scaffold_' in line or 'super_' in line or 'ta_' in line:
                    scaffoldOutList += [list(outputTuple)] # separate outputs for scaffolds
                else:
                    if chrCount >= 25: # if chromosome count above 24, use random coloring
                        outputFile.write('chr - %s %s 0 %s %d,%d,%d\n' % tuple(list(outputTuple)[:3]+
                                                                [randint(1,255),randint(1,255),randint(1,255)]))
                    else: # if chromosome count below 25, use predictable coloring by chromosome number
                        outputFile.write('chr - %s %s 0 %s chr%d\n' % outputTuple)
                        chrCount += 1 # add one to chromosome count
            if chrCount == 25: #FIXME ADD MORE THAN 25 CHROMOSOMES; this is fixed I hope, will output as warning
                print('Warning: Chromosome count exceeded in %s (anymore chromosomes may not generate useful display...'
                      %outputFilename)
                break
        # WORKING ON THIS!!!
        # if scaffolds/super/etc exist
        if scaffoldOutList:
            try:
                # sort the scaffolds in order by number
                scaffoldOutList.sort(key = lambda x: int(x[0].split('_')[-1]))
            except:
                # try running analysis without the file if spits out an error, if that doesn't work, contact Joshua Levy
                print 'There is an error in '+ file
                exit()
                # output each scaffold to karyotype file
            for scaffold in scaffoldOutList:
                scaffold[-1] = chrCount # chrCount used to color scaffolds in circos plot
                if 'ta' in scaffold[0]: # wheat has generated problems in past...
                    outputFile.write('chr - %s %s 0 %s %d,%d,%d\n' % tuple(scaffold[:3]+[randint(1,255),randint(1,255),
                                                                                      randint(1,255)]))
                # if chromosomes less than 25, write to file
                elif chrCount < 25:
                    outputFile.write('chr - %s %s 0 %s chr%d\n' % tuple(scaffold))
                    chrCount += 1
                else: # if chromosomes number greater than 24 and if not query, then randomize color
                    outputFile.write('chr - %s %s 0 %s %d,%d,%d\n' % tuple(scaffold[:3]+[randint(1,255),
                                                                            randint(1,255),randint(1,255)]))
        # close fai and karyotype files
        outputFile.close()
        inputFile.close()





# return a list of the karyotype files
    return outputFilename


    # edit file to liking afterwards...


    """Bd1	75071545	5	80	81
    Bd1_centromere_containing_Bradi1g41430	46184	76009985	80	81
    Bd2	59130575	76056752	80	81
    Bd3	59640145	135926465	80	81
    Bd4	48594894	196312117	80	81
    Bd5	28630136	245514453	80	81
    scaffold_12	23566	274502479	80	81
    scaffold_135	3881	274526354	80	81
    scaffold_14	20560	274530297	80	81
    scaffold_180	1933	274551128	80	81"""

    """chr - hsX X 0 153692391 chrx
    chr - hsY Y 0 50286555 chry
    chr - hs1 1 0 246127941 chr1
    chr - hs2 2 0 243615958 chr2
    chr - hs3 3 0 199344050 chr3
    chr - hs4 4 0 191731959 chr4
    chr - hs5 5 0 181034922 chr5
    chr - hs6 6 0 170914576 chr6
    chr - hs7 7 0 158545518 chr7
    chr - hs8 8 0 146308819 chr8
    chr - hs9 9 0 136372045 chr9
    chr - hs10 10 0 135037215 chr10
    chr - hs11 11 0 134482954 chr11
    chr - hs12 12 0 132078379 chr12
    chr - hs13 13 0 113042980 chr13
    chr - hs14 14 0 105311216 chr14
    chr - hs15 15 0 100256656 chr15
    chr - hs16 16 0 90041932 chr16"""

    """Bstacei_316_v1.0.softmasked.fa.fai
    Osativa_323_v7.0.softmasked.fa.fai
    Phallii_308_v2.0.softmasked.fa.fai
    Pvirgatum_383_v3.0.softmasked.fa.fai
    Sbicolor_313_v3.0.softmasked.fa.fai
    Sitalica_312_v2.softmasked.fa.fai
    Sviridis_311_v1.0.softmasked.fa.fai"""