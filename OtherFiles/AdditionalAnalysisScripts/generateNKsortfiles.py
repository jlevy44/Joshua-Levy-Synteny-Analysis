sortFile = 't.PAC2_0.383.sort'

open('t.PAC2_0.383K.sort','w').close()
open('t.PAC2_0.383N.sort','w').close()


sortFileInput = open('t.PAC2_0.383.sort','r')
sortFileK = open('t.PAC2_0.383K.sort','w')
sortFileN = open('t.PAC2_0.383N.sort','w')
for line in sortFileInput:
    if 'K' in line:
        sortFileK.write(line)
    if 'N' in line:
        a=1

#if 'J' in line:
