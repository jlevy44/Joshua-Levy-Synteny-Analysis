import os, sys,subprocess
"""
grep "BhD" q.PAC4GC.564.gff3  > transformed_GFFFiles/q.PAC4GC.001.gff3

grep "BhS" q.PAC4GC.564.gff3   > transformed_GFFFiles/q.PAC4GC.002.gff3

grep "BhD" q.PAC4GC.590.gff3   > transformed_GFFFiles/q.PAC4GC.003.gff3

grep "BhS" q.PAC4GC.590.gff3  > transformed_GFFFiles/q.PAC4GC.004.gff3

grep "BhD" q.PAC4GC.601.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.005.gff3

grep "BhS" q.PAC4GC.601.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.006.gff3

grep "BhD" q.PAC4GC.598.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.007.gff3

grep "BhS" q.PAC4GC.598.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.008.gff3

grep "BhD" q.PAC4GC.596.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.009.gff3

grep "BhS" q.PAC4GC.596.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.010.gff3

grep "BhD" q.PAC4GC.593.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.011.gff3

grep "BhS" q.PAC4GC.593.gff3 | sed 's/ABR113./ABR113/g' > transformed_GFFFiles/q.PAC4GC.012.gff3
"""
"""change fasta files"""

"""
sed 's/Bh1D/BhD/g' BhybridumABR113D_001_v1.0.softmasked.fa > BhybridumABR113D_001_v1.0.softmasked.fa
sed 's/Bh1D/BhD/g' BhybridumABR113D_001_v1.0.softmasked.fa.fai > BhybridumABR113D_001_v1.0.softmasked.fa.fai
sed 's/Bh2D/BhD/g' BhybridumBd28D_003_v1.0.softmasked.fa > BhybridumBd28D_003_v1.0.softmasked.fa
sed 's/Bh2D/BhD/g' BhybridumBd28D_003_v1.0.softmasked.fa.fai > BhybridumBd28D_003_v1.0.softmasked.fa.fai
sed 's/Bh3D/BhD/g' BhybridumBdTR6gD_005_v1.0.softmasked.fa > BhybridumBdTR6gD_005_v1.0.softmasked.fa
sed 's/Bh3D/BhD/g' BhybridumBdTR6gD_005_v1.0.softmasked.fa.fai > BhybridumBdTR6gD_005_v1.0.softmasked.fa.fai
sed 's/Bh6D/BhD/g' BhybridumBhyb30D_011_v1.0.softmasked.fa > BhybridumBhyb30D_011_v1.0.softmasked.fa
sed 's/Bh6D/BhD/g' BhybridumBhyb30D_011_v1.0.softmasked.fa.fai > BhybridumBhyb30D_011_v1.0.softmasked.fa.fai
sed 's/Bh5D/BhD/g' BhybridumBhyb36D_009_v1.0.softmasked.fa > BhybridumBhyb36D_009_v1.0.softmasked.fa
sed 's/Bh5D/BhD/g' BhybridumBhyb36D_009_v1.0.softmasked.fa.fai > BhybridumBhyb36D_009_v1.0.softmasked.fa.fai
sed 's/Bh4D/BhD/g' BhybridumPob1D_007_v1.0.softmasked.fa > BhybridumPob1D_007_v1.0.softmasked.fa
sed 's/Bh4D/BhD/g' BhybridumPob1D_007_v1.0.softmasked.fa.fai > BhybridumPob1D_007_v1.0.softmasked.fa.fai
"""



listUnOut = str(subprocess.Popen(['ls','.'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        .stdout.read()).split('\n')
speciesDictionary = {'590':['003','004'],'601':['005','006'],'598':['007','008'],'596':['009','010'],'593':['011','012']}
for unoutfile in listUnOut:
    if unoutfile and unoutfile.endswith('.unout'):
        unoutTarget,unoutQuery = tuple(unoutfile.split('-'))
        if 'BhD' in unoutTarget:
            unoutTarget = unoutTarget.replace('564','001').replace('.pre.BhD','')
        else:
            unoutTarget = unoutTarget.replace('564','002').replace('.pre.BhS','')
        if 'BhD' in unoutQuery:
            unoutQuery = unoutQuery.replace(unoutQuery[unoutQuery.find('.')+1:unoutQuery.find('.pre')],
                                            speciesDictionary[unoutQuery[unoutQuery.find('.')+1:unoutQuery.find('.pre')]][0]).replace('.pre.BhD','')
        elif 'BhS' in unoutQuery:
            unoutQuery = unoutQuery.replace(unoutQuery[unoutQuery.find('.') + 1:unoutQuery.find('.pre')],
                                            speciesDictionary[unoutQuery[unoutQuery.find('.') + 1:unoutQuery.find('.pre')]][1]).replace('.pre.BhS', '')
        input = open(unoutfile,'r')
        outputText = input.read().replace('ABR113.','ABR113')
        outputFile = '-'.join(unout for unout in [unoutTarget,unoutQuery])
        open(outputFile, 'w').close()
        outputF = open(outputFile,'w')
        outputF.write(outputText)
        outputF.close()
        input.close()
