import subprocess

circosConfigFilesPath = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/' \
                        'CircosInputs/ConfigFiles/'
circosOutPath = '/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/DNA_Segment_Extraction/' \
                'CircosInputs/Data'

speciesList = ['523', '524-Hist']
subprocess.call(['/Applications/circos-0.69/bin/circos', '-conf', circosConfigFilesPath + 'circos.conf', '-outputfile',
                 '%s-%s' % (speciesList[0], speciesList[1]), '-outputdir',circosOutPath])