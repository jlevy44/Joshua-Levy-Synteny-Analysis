import os, sys
os.chdir('v0/Bdist_%s_v0'%(sys.argv[1]))
for file in os.listdir('.'):
    if file.startswith('data'):
        os.remove(file)
try:
    for sample in sys.argv[2:]:
        os.chdir('../Bdist_%s_v0'%(sample))
        for file in os.listdir('.'):
            if file.startswith('data') or file.startswith('Odata'):
                os.remove(file)
except:
    pass