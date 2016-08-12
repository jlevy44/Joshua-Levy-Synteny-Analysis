import subprocess, os, shutil

subprocess.call('pwd')
#os.chdir('allMapsJoshTest')
try:
    listFiles=str(subprocess.Popen(['ls','allMapsJoshTest'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.read()).split('\n')
    for file in listFiles:
        if file.endswith('.bed'):
            shutil.move('/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/ALLMAPS/allMapsJoshTest/'+file,'/Users/jlevy/Documents/Joshua_Research_JGI/Joshua_Levy_Research_Outputs/ALLMAPS/allMapsJoshTest/BdFiles')
except:
    print 'Error: Unable to move output files. Move *.pdf to ALLMAPS_Chromosome_Images and make sure BdFiles'
#shutil.move('*.bed','./allMapsJoshTest/BdFiles')
#subprocess.call(['python','createNewGenome.py'])