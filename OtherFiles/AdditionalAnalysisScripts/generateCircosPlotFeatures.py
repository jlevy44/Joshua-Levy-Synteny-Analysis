import subprocess

outputPlotFilename = '283.323.histogram.'

open(outputPlotFilename+'ShellRatio.txt','w').close()
open(outputPlotFilename+'NonSyntenicRatio.txt','w').close()

shellFileOut = open(outputPlotFilename+'ShellRatio.txt','w')
nonSynFileOut = open(outputPlotFilename+'NonSyntenicRatio.txt','w')
inputFile = open('all.shell.core.TE.low_depth.SNPs.recombRate.TEmatrixINS.TEmatrixABS.synteny.filt.bedgraph.ShellRatio'
                 '.tr.txt','r')
read = 0
for line in inputFile:
    if line:
        if read == 1:
            lineList = line.split('\t')
            shellFileOut.write('%s %s %s %s\n' %(lineList[0],lineList[1],lineList[2],lineList[-2]))
            nonSynFileOut.write('%s %s %s %s\n' %(lineList[0],lineList[1],lineList[2],lineList[-1]))


        if 'start' in line:
            read = 1

    else:
        break

inputFile.close()
shellFileOut.close()
nonSynFileOut.close()
circosConfigFilesPath = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/' \
                        'BdistachyonAnalysis/CircosInputs/ConfigFiles/'
circosOutPath = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/' \
                'BdistachyonAnalysis/CircosInputs/Data/Plots'

speciesList = ['283','323']
subprocess.call(['/Applications/circos-0.69/bin/circos', '-conf', circosConfigFilesPath + 'circos.conf', '-outputfile',
                 '%s-%s-%s' % (speciesList[0], speciesList[1],'WithShellRatio-NonSyntenic-Plots'), '-outputdir', circosOutPath])