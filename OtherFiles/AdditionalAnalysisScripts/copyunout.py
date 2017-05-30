import shutil

unouts = """PAC4GC.001-PAC2_0.275_5.unout
PAC4GC.001-PAC2_0.285_5.unout
PAC4GC.001-PAC2_0.308_5.unout
PAC4GC.001-PAC2_0.311_5.unout
PAC4GC.001-PAC2_0.312_5.unout
PAC4GC.001-PAC2_0.313_5.unout
PAC4GC.001-PAC2_0.314_5.unout
PAC4GC.001-PAC2_0.316_5.unout
PAC4GC.001-PAC2_0.323_5.unout
PAC4GC.001-PAC2_0.348_5.unout
PAC4GC.001-PAC2_0.354_5.unout
PAC4GC.001-PAC2_0.362_5.unout
PAC4GC.001-PAC2_0.383_5.unout
PAC4GC.001-PAC4GC.003_5.unout
PAC4GC.001-PAC4GC.004_5.unout
PAC4GC.001-PAC4GC.005_5.unout
PAC4GC.001-PAC4GC.006_5.unout
PAC4GC.001-PAC4GC.007_5.unout
PAC4GC.001-PAC4GC.008_5.unout
PAC4GC.001-PAC4GC.009_5.unout
PAC4GC.001-PAC4GC.010_5.unout
PAC4GC.001-PAC4GC.011_5.unout
PAC4GC.001-PAC4GC.012_5.unout"""
allowed = ['003','005','007','009','011','348','314','362','354','316','323']
unoutList = unouts.split()
for unout in unoutList:
    if unout:
        for species in allowed:
            if species in unout:
                shutil.copy(unout,'/global/projectb/scratch/jlevy/Spring2017/Synteny/Bhybridum/D/UnOutFiles/')