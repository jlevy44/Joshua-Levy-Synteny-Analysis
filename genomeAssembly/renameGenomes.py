import sys, shutil, os
from collections import defaultdict

e = sys.argv[1:]
counter = 1
fastaCorrespond = defaultdict(list)
os.mkdir('v0')
for file in os.listdir('.'):
    if file.endswith('.fa'):
        os.mkdir('v0/%s_%s_v0'%(e[0],'0'* (3-len(str(counter)))+ str(counter)))
        shutil.move(file,'v0/%s_%s_v0/%s_%s_v0.fa'%(e[0],'0'* (3-len(str(counter)))+ str(counter),e[0],'0'* (3-len(str(counter)))+ str(counter)))
        fastaCorrespond[file] = '%s_%s_v0.fa'%(e[0],'0'* (3-len(str(counter)))+ str(counter))
        counter += 1

with open('CorrespondenceAssemblies.txt','w') as f:
    f.write('\n'.join(original+'\t'+fastaCorrespond[original] for original in fastaCorrespond.keys()))

