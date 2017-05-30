import shutil

convertText = """AmTr_v1.0_scaffolds.hardmasked2gff.gff 291
Bd18-1.h.repeatMasked.gff 362
Bd28.genome.DS.allmaps.chrom.repeatMasked.noSubHeaders.BhD.gff 003
Bd28.genome.DS.allmaps.chrom.repeatMasked.noSubHeaders.BhS.gff 004
BdTR3C.h.repeatMasked.gff 354
BdTR6g_PNWH.DS.chr.mapMer.BhU.repeatMasked.noSubHeaders.BhD.gff 005
BdTR6g_PNWH.DS.chr.mapMer.BhU.repeatMasked.noSubHeaders.BhS.gff 006
BdTR8i.h.repeatMasked.gff 348
Bhyb30_AOZWG.mainGenome.BhU.repeatMasked.noSubHeaders.BhD.gff 011
Bhyb30_AOZWG.mainGenome.BhU.repeatMasked.noSubHeaders.BhS.gff 012
Bhyb36_AOZWH.mainGenome.BhU.repeatMasked.noSubHeaders.BhD.gff 009
Bhyb36_AOZWH.mainGenome.BhU.repeatMasked.noSubHeaders.BhS.gff 010
Brachypodium_distachyon.mainGenome.repeatMasked.gff 314
Brachypodium_hybridum.mainGenome.scaffolds.gapfilled.091816.repeatMasked.noSubHeaders.BhD.gff 001
Brachypodium_hybridum.mainGenome.scaffolds.gapfilled.091816.repeatMasked.noSubHeaders.BhS.gff 002
Brachypodium_stacei_ABR114.mainGenome.repeatMasked.gff 316
Othomaeum_386_v1.0.repeatMasked.gff 386
Pob1.DS.chr.mapMer.BhU.repeatMasked.noSubHeaders.BhD.gff 006
Pob1.DS.chr.mapMer.BhU.repeatMasked.noSubHeaders.BhS.gff 007
Sorghum_bicolor.main_genome.scaffolds.repeatMasked.gff 313
Sp_repeats.gff3 290
Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.repeats.from.softmasked.genome.gff 443
Zmarina_324_v2.2.repeatMasked.gff 324
maize.repeat_region.transformedFrom.bed.gff2 284
musa.RepeatMasker.fromEnsemblgff.gff3 304
pineapple.20150427.xMasked.gff2 321
rice_osa1r7_rm.gff3 323
"""

for line in convertText.split('\n'):
    if line:
        infile = line.split()[0]
        outfile = line.split()[1].strip('\n')+'_transposonDensity'+infile[infile.rfind('.'):]
        shutil.move(infile,outfile)