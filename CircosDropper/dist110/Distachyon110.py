import re
from ete2 import Tree
import pandas as pd
import numpy as np
oldDist = {line.split('_')[0].replace('Bdistachyon','').lower() : line.split('_')[1] for line in """BdistachyonABR2_337_v1.softmasked.fa.fai
BdistachyonABR3_343_v1.softmasked.fa.fai
BdistachyonABR4_364_v1.softmasked.fa.fai
BdistachyonABR5_379_v1.softmasked.fa.fai
BdistachyonABR6_336_v1.softmasked.fa.fai
BdistachyonABR7_369_v1.softmasked.fa.fai
BdistachyonABR8_356_v1.softmasked.fa.fai
BdistachyonABR9_333_v1.softmasked.fa.fai
BdistachyonAdi2_359_v1.softmasked.fa.fai
BdistachyonAdi_372_v1.softmasked.fa.fai
BdistachyonArn1_355_v1.softmasked.fa.fai
BdistachyonBd18_362_v1.softmasked.fa.fai
BdistachyonBd1_349_v1.softmasked.fa.fai
BdistachyonBd21AnntCtrl_331_v2.0.softmasked.fa.fai
BdistachyonBd21AsmbCtrl_361_v1.softmasked.fa.fai
BdistachyonBd21_378_v1.softmasked.fa.fai
BdistachyonBd21_460_v1.0.softmasked.fa.fai
BdistachyonBd21v283_v2.0.softmasked.fa.fai
BdistachyonBd29_346_v1.softmasked.fa.fai
BdistachyonBd2_353_v1.softmasked.fa.fai
BdistachyonBd30_344_v1.softmasked.fa.fai
BdistachyonBdTR10c_374_v1.softmasked.fa.fai
BdistachyonBdTR11a_380_v1.softmasked.fa.fai
BdistachyonBdTR11g_357_v1.softmasked.fa.fai
BdistachyonBdTR11i_363_v1.softmasked.fa.fai
BdistachyonBdTR12c_352_v1.softmasked.fa.fai
BdistachyonBdTR13a_334_v1.softmasked.fa.fai
BdistachyonBdTR13c_365_v1.softmasked.fa.fai
BdistachyonBdTR1i_345_v1.softmasked.fa.fai
BdistachyonBdTR2b_376_v1.softmasked.fa.fai
BdistachyonBdTR2g_367_v1.softmasked.fa.fai
BdistachyonBdTR3c_354_v1.softmasked.fa.fai
BdistachyonBdTR5i_370_v1.softmasked.fa.fai
BdistachyonBdTR7a_329_v1.softmasked.fa.fai
BdistachyonBdTR8i_348_v1.softmasked.fa.fai
BdistachyonBdTR9k_358_v1.softmasked.fa.fai
BdistachyonBis_338_v1.softmasked.fa.fai
BdistachyonFoz1_366_v1.softmasked.fa.fai
BdistachyonJer1_375_v1.softmasked.fa.fai
BdistachyonKah_342_v1.softmasked.fa.fai
BdistachyonKoz_330_v1.softmasked.fa.fai
BdistachyonKoz_368_v1.softmasked.fa.fai
BdistachyonLuc1_347_v1.softmasked.fa.fai
BdistachyonMig3_377_v1.softmasked.fa.fai
BdistachyonMon3_350_v1.softmasked.fa.fai
BdistachyonMur1_373_v1.softmasked.fa.fai
BdistachyonPer1_326_v1.softmasked.fa.fai
BdistachyonRon2_360_v1.softmasked.fa.fai
BdistachyonS8iiC_340_v1.softmasked.fa.fai
BdistachyonSig2_371_v1.softmasked.fa.fai
BdistachyonTek_341_v1.softmasked.fa.fai
BdistachyonUni2_327_v1.softmasked.fa.fai
Bdistachyon_314_v3.0.softmasked.fa.fai""".splitlines() }

more = { line.split()[0].lower(): line.split()[1] for line in """Koz-1 330
Tek-2 341
ABR6_r 336
ABR9_r 333
Bd2-3 353
Adi-12 359
Bd18-1 362
Koz-3 368
Bd21-3_r 460
Bd30-1 344
Bd1-1 349
Kah-5 342
Bis-1 338
Bd21Control 361
Adi-2 372
""".splitlines()}

bad = list(set("""356
381
328
339
351
335
332""".splitlines() + [line.split()[1] for line in """9     681
11    374
13    380
16    346
36    646
49    673
75    370
81    367
91    708
95    669
99    358""".splitlines()]))

list110 = list(set("""283
326
327
329
330
331
333
334
336
337
338
340
341
342
343
344
345
346
347
348
349
350
352
353
354
355
356
357
358
359
360
361
362
363
364
365
366
367
368
369
370
371
372
373
374
375
376
377
379
380
460
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666
667
668
669
670
671
672
673
674
675
676
677
678
679
680
681
682
683
684
685
686
687
688
689
690
691
692
693
694
695
696
697
698
699
700
701
702
703
704
705
706
707
708
314
""".splitlines())-set(bad))
inputStr = """646     Brachypodium distachyon 018     lo.v1.1 31081
647     Brachypodium distachyon 033     lo.v1.1 30664
648     Brachypodium distachyon 023     lo.v1.1 31298
649     Brachypodium distachyon 053     lo.v1.1 31664
650     Brachypodium distachyon 035     lo.v1.1 31706
651     Brachypodium distachyon 010     lo.v1.1 31596
652     Brachypodium distachyon 050     lo.v1.1 30326
653     Brachypodium distachyon 044     lo.v1.1 31817
654     Brachypodium distachyon 056     lo.v1.1 31375
655     Brachypodium distachyon 060     lo.v1.1 30321
656     Brachypodium distachyon 006     lo.v1.1 29422
657     Brachypodium distachyon 014     lo.v1.1 31415
658     Brachypodium distachyon 005     lo.v1.1 31263
659     Brachypodium distachyon 017     lo.v1.1 31460
660     Brachypodium distachyon 007     lo.v1.1 31260
661     Brachypodium distachyon 011     lo.v1.1 31192
662     Brachypodium distachyon 036     lo.v1.1 31415
663     Brachypodium distachyon 012     lo.v1.1 31673
664     Brachypodium distachyon 052     lo.v1.1 31236
665     Brachypodium distachyon 040     lo.v1.1 31833
666     Brachypodium distachyon 032     lo.v1.1 31546
667     Brachypodium distachyon 022     lo.v1.1 31205
668     Brachypodium distachyon 003     lo.v1.1 32319
669     Brachypodium distachyon 046     lo.v1.1 27790
670     Brachypodium distachyon 062     lo.v1.1 30346
671     Brachypodium distachyon 026     lo.v1.1 29986
672     Brachypodium distachyon 037     lo.v1.1 30374
673     Brachypodium distachyon 059     lo.v1.1 29698
674     Brachypodium distachyon 015     lo.v1.1 30942
675     Brachypodium distachyon 019     lo.v1.1 31245
676     Brachypodium distachyon 055     lo.v1.1 30651
677     Brachypodium distachyon 063     lo.v1.1 31396
678     Brachypodium distachyon 049     lo.v1.1 31153
679     Brachypodium distachyon 001     lo.v1.1 29602
680     Brachypodium distachyon 041     lo.v1.1 30652
681     Brachypodium distachyon 047     lo.v1.1 28148
682     Brachypodium distachyon 034     lo.v1.1 31282
683     Brachypodium distachyon 028     lo.v1.1 28892
684     Brachypodium distachyon 057     lo.v1.1 30715
685     Brachypodium distachyon 013     lo.v1.1 30239
686     Brachypodium distachyon 042     lo.v1.1 29757
687     Brachypodium distachyon 058     lo.v1.1 30484
688     Brachypodium distachyon 054     lo.v1.1 30757
689     Brachypodium distachyon 024     lo.v1.1 29692
690     Brachypodium distachyon 045     lo.v1.1 29977
691     Brachypodium distachyon 048     lo.v1.1 29741
692     Brachypodium distachyon 025     lo.v1.1 29880
693     Brachypodium distachyon 008     lo.v1.1 31285
694     Brachypodium distachyon 030     lo.v1.1 30598
695     Brachypodium distachyon 031     lo.v1.1 30775
696     Brachypodium distachyon 002     lo.v1.1 30275
697     Brachypodium distachyon 051     lo.v1.1 30858
698     Brachypodium distachyon 064     lo.v1.1 30872
699     Brachypodium distachyon 016     lo.v1.1 31445
700     Brachypodium distachyon 009     lo.v1.1 30872
701     Brachypodium distachyon 020     lo.v1.1 31998
702     Brachypodium distachyon 004     lo.v1.1 30958
703     Brachypodium distachyon 039     lo.v1.1 31330
704     Brachypodium distachyon 029     lo.v1.1 31101
705     Brachypodium distachyon 021     lo.v1.1 31465
706     Brachypodium distachyon 038     lo.v1.1 30706
707     Brachypodium distachyon 027     lo.v1.1 31253
708     Brachypodium distachyon 061     lo.v1.1 24880"""
traceback = {line.split()[3].lower():line.split()[0] for line in inputStr.splitlines()}
tree_raw = """(100_0000007[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,BdtwD[&Created="Tue May 23 17:01:31 CEST 2017"]:0.001437,(((Ibd-189[&Created="Tue May 23 17:01:31 CEST 2017"]:0.034563,(Ibd-107[&Created="Tue May 23 17:01:31 CEST 2017"]:0.05466,Bd28[&Created="Tue May 23 17:01:31 CEST 2017"]:0.036353)[&"Consensus support(%)"=100.0]:0.049391)[&"Consensus support(%)"=79.0]:0.030434,((BdTR6g[&Created="Tue May 23 17:01:31 CEST 2017"]:9.8E-5,BdTRgD[&Created="Tue May 23 17:01:31 CEST 2017"]:0.002115)[&"Consensus support(%)"=100.0]:0.056387,((ABR113[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,ABRD[&Created="Tue May 23 17:01:31 CEST 2017"]:1.2E-5)[&"Consensus support(%)"=100.0]:0.092386,((Bhyb35[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,BhybzeroD[&Created="Tue May 23 17:01:31 CEST 2017"]:0.001528)[&"Consensus support(%)"=100.0]:0.111374,((Pob1[&Created="Tue May 23 17:01:31 CEST 2017"]:1.3E-5,PobD[&Created="Tue May 23 17:01:31 CEST 2017"]:0.002055)[&"Consensus support(%)"=100.0]:0.103278,(('Ita-Sic-LPA3.2'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.021188,(Ita-Sic-CSR-7[&Created="Tue May 23 17:01:31 CEST 2017"]:0.003455,Ita-Sic-CSR-6[&Created="Tue May 23 17:01:31 CEST 2017"]:0.009187)[&"Consensus support(%)"=100.0]:0.0114)[&"Consensus support(%)"=100.0]:0.078601,(('Spa-Sou-GR6.4'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.039689,(Arn1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.023585,(Spa-Nor-S6D-5[&Created="Tue May 23 17:01:31 CEST 2017"]:5.7E-5,Mon3[&Created="Tue May 23 17:01:31 CEST 2017"]:9.0E-6)[&"Consensus support(%)"=100.0]:0.045961)[&"Consensus support(%)"=100.0]:0.024575)[&"Consensus support(%)"=100.0]:0.085391,((((Bd30-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.064361,(('Spa-Sou-J6.2'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.060289,'Spa-Sou-AB1.4-2'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.064205)[&"Consensus support(%)"=97.0]:0.004785,('Spa-SouGR3.6-5'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.052727,'Spa-Sou-J4.3.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.063768)[&"Consensus support(%)"=97.0]:0.004819)[&"Consensus support(%)"=77.0]:0.003156)[&"Consensus support(%)"=72.0]:0.005786,((ABR2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.07839,((Luc1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.068398,(Spa-Nor-S19A[&Created="Tue May 23 17:01:31 CEST 2017"]:0.068487,((RON2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.060375,(Per1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.051742,ABR6_r[&Created="Tue May 23 17:01:31 CEST 2017"]:0.053389)[&"Consensus support(%)"=100.0]:0.011299)[&"Consensus support(%)"=97.0]:0.006664,((ABR3[&Created="Tue May 23 17:01:31 CEST 2017"]:0.062893,ABR5[&Created="Tue May 23 17:01:31 CEST 2017"]:0.056695)[&"Consensus support(%)"=88.0]:0.00424,(Mur1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.062571,(Foz1[&Created="Tue May 23 17:01:31 CEST 2017"]:8.0E-6,Sig2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=100.0]:0.062345)[&"Consensus support(%)"=77.0]:0.00291)[&"Consensus support(%)"=56.00000000000001]:0.002024)[&"Consensus support(%)"=35.0]:0.001473)[&"Consensus support(%)"=65.0]:0.002631)[&"Consensus support(%)"=77.0]:0.003115,(((ABR4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.07614,(Spa-Nor-S11B-2[&Created="Tue May 23 17:01:31 CEST 2017"]:3.5E-5,'Spa-Nor-S11A.3'[&Created="Tue May 23 17:01:31 CEST 2017"]:3.0E-6)[&"Consensus support(%)"=100.0]:0.07481)[&"Consensus support(%)"=62.0]:0.003257,((Spa-Nor-S12B-4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.043396,Mig3[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=100.0]:0.016012,(Jer1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.053535,('Spa-Sou-HU3.4'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.061167,('Spa-Nor-S22C.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,S8iiC[&Created="Tue May 23 17:01:31 CEST 2017"]:1.4E-5)[&"Consensus support(%)"=100.0]:0.057868)[&"Consensus support(%)"=50.0]:0.003081)[&"Consensus support(%)"=92.0]:0.004435)[&"Consensus support(%)"=22.0]:0.001384)[&"Consensus support(%)"=34.0]:0.002433,((Spa-Nor-S6E[&Created="Tue May 23 17:01:31 CEST 2017"]:2.0E-6,Spa-Nor-S6B-3[&Created="Tue May 23 17:01:31 CEST 2017"]:1.15E-4)[&"Consensus support(%)"=100.0]:0.062232,(Uni2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.013575,((Spa-Nor-S11D-2[&Created="Tue May 23 17:01:31 CEST 2017"]:5.0E-6,Spa-Nor-S17D-1[&Created="Tue May 23 17:01:31 CEST 2017"]:7.2E-5)[&"Consensus support(%)"=100.0]:0.060787,(Spa-Nor-S16D[&Created="Tue May 23 17:01:31 CEST 2017"]:0.06069,Spa-Nor-S22B-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.054807)[&"Consensus support(%)"=56.99999999999999]:0.003096)[&"Consensus support(%)"=66.0]:0.002635)[&"Consensus support(%)"=56.00000000000001]:0.002823)[&"Consensus support(%)"=31.0]:0.002151)[&"Consensus support(%)"=63.0]:0.00199)[&"Consensus support(%)"=100.0]:0.008509)[&"Consensus support(%)"=36.0]:0.002014,((ABR7[&Created="Tue May 23 17:01:31 CEST 2017"]:0.040409,('Spa-Sou-Sg2.1.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.062679,('Spa-Sou-AB3.3'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.066673,'Spa-Sou-CU1.6-1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.066007)[&"Consensus support(%)"=84.0]:0.004113)[&"Consensus support(%)"=52.0]:0.005232)[&"Consensus support(%)"=19.0]:0.002081,('Spa-Sou-CU5.5'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.059488,('Spa-Sou-Z3.6'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.059035,Spa-GU5-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.07555)[&"Consensus support(%)"=56.99999999999999]:0.003181)[&"Consensus support(%)"=12.0]:0.001511)[&"Consensus support(%)"=28.000000000000004]:0.002972)[&"Consensus support(%)"=49.0]:0.007456)[&"Consensus support(%)"=100.0]:0.076201,((Bd18-1[&Created="Tue May 23 17:01:31 CEST 2017"]:1.7E-5,BdEighteen[&Created="Tue May 23 17:01:31 CEST 2017"]:6.0E-5)[&"Consensus support(%)"=100.0]:0.087862,((Arm-Arm3i[&Created="Tue May 23 17:01:31 CEST 2017"]:0.047191,(Arm-Arm2B[&Created="Tue May 23 17:01:31 CEST 2017"]:0.057943,(((Arm-Arm3G[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,Arm-Arm3A[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=100.0]:0.033892,(Geo-G34i6[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,Geo-G34i2[&Created="Tue May 23 17:01:31 CEST 2017"]:1.16E-4)[&"Consensus support(%)"=100.0]:0.026945)[&"Consensus support(%)"=100.0]:0.025688,(Geo-G32i2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.040857,Geo-G31i4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.034648)[&"Consensus support(%)"=100.0]:0.009976)[&"Consensus support(%)"=100.0]:0.01758)[&"Consensus support(%)"=100.0]:0.011045)[&"Consensus support(%)"=100.0]:0.038414,(((Gaz-7[&Created="Tue May 23 17:01:31 CEST 2017"]:0.014983,Gaz-2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.006309)[&"Consensus support(%)"=100.0]:0.048483,(Gaz-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.040727,(Gaz-8[&Created="Tue May 23 17:01:31 CEST 2017"]:0.031026,(BdTR1i[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR2G[&Created="Tue May 23 17:01:31 CEST 2017"]:8.5E-5,BdTR2B[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=94.0]:4.6E-5)[&"Consensus support(%)"=100.0]:0.028769)[&"Consensus support(%)"=100.0]:0.040994)[&"Consensus support(%)"=100.0]:0.014306)[&"Consensus support(%)"=100.0]:0.022901,(Koz-3[&Created="Tue May 23 17:01:31 CEST 2017"]:0.073524,((Adi-10[&Created="Tue May 23 17:01:31 CEST 2017"]:0.09382,((Bd3-1_r[&Created="Tue May 23 17:01:31 CEST 2017"]:0.07902,(Bd2-3[&Created="Tue May 23 17:01:31 CEST 2017"]:0.089644,(Bd21Control[&Created="Tue May 23 17:01:31 CEST 2017"]:0.121294,Bdistachyon[&Created="Tue May 23 17:01:31 CEST 2017"]:0.121294,(Vasc1[&Created="Tue May 23 17:01:31 CEST 2017"]:5.0E-6,Bd21-3_r[&Created="Tue May 23 17:01:31 CEST 2017"]:0.001564)[&"Consensus support(%)"=100.0]:0.092432)[&"Consensus support(%)"=100.0]:0.015726)[&"Consensus support(%)"=100.0]:0.022514)[&"Consensus support(%)"=100.0]:0.030548,(Adi-9[&Created="Tue May 23 17:01:31 CEST 2017"]:0.07603,(((Kah-6[&Created="Tue May 23 17:01:31 CEST 2017"]:0.082505,Kah-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.063036)[&"Consensus support(%)"=94.0]:0.014497,(Kah-5[&Created="Tue May 23 17:01:31 CEST 2017"]:0.073166,(BdTR5A[&Created="Tue May 23 17:01:31 CEST 2017"]:5.0E-6,BdTR5I[&Created="Tue May 23 17:01:31 CEST 2017"]:1.48E-4)[&"Consensus support(%)"=100.0]:0.076136)[&"Consensus support(%)"=71.0]:0.007694)[&"Consensus support(%)"=72.0]:0.005691,((Adi-15[&Created="Tue May 23 17:01:31 CEST 2017"]:0.04534,(BdTR11E[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR11G[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR11I[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,BdTR11A[&Created="Tue May 23 17:01:31 CEST 2017"]:7.6E-5)[&"Consensus support(%)"=56.00000000000001]:4.9E-5)[&"Consensus support(%)"=56.99999999999999]:3.8E-5)[&"Consensus support(%)"=100.0]:0.03159)[&"Consensus support(%)"=100.0]:0.042963,((((Adi-4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.050504,Adi-2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.054959)[&"Consensus support(%)"=100.0]:0.007771,(Adi-12[&Created="Tue May 23 17:01:31 CEST 2017"]:0.024141,(BdTR9M[&Created="Tue May 23 17:01:31 CEST 2017"]:1.4E-5,BdTR9K[&Created="Tue May 23 17:01:31 CEST 2017"]:1.47E-4)[&"Consensus support(%)"=100.0]:0.027932)[&"Consensus support(%)"=100.0]:0.02992)[&"Consensus support(%)"=100.0]:0.014363,(BdTR12B[&Created="Tue May 23 17:01:31 CEST 2017"]:3.2E-5,BdTR12c[&Created="Tue May 23 17:01:31 CEST 2017"]:1.0E-5)[&"Consensus support(%)"=100.0]:0.048067)[&"Consensus support(%)"=100.0]:0.033272,(BdTR10D[&Created="Tue May 23 17:01:31 CEST 2017"]:1.85E-4,BdTR10C[&Created="Tue May 23 17:01:31 CEST 2017"]:1.5E-5)[&"Consensus support(%)"=100.0]:0.08515)[&"Consensus support(%)"=54.0]:0.005083)[&"Consensus support(%)"=38.0]:0.003783)[&"Consensus support(%)"=44.0]:0.003987)[&"Consensus support(%)"=93.0]:0.006615)[&"Consensus support(%)"=70.0]:0.003962)[&"Consensus support(%)"=98.0]:0.011711,((BdTR13B[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR13N[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR13a[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,(BdTR13C[&Created="Tue May 23 17:01:31 CEST 2017"]:1.16E-4,Bis-1[&Created="Tue May 23 17:01:31 CEST 2017"]:3.0E-6)[&"Consensus support(%)"=52.0]:1.26E-4)[&"Consensus support(%)"=42.0]:7.0E-5)[&"Consensus support(%)"=42.0]:3.6E-5)[&"Consensus support(%)"=100.0]:0.069962,(Koz-5[&Created="Tue May 23 17:01:31 CEST 2017"]:0.06546,((BdTRc[&Created="Tue May 23 17:01:31 CEST 2017"]:3.8E-4,(BdTR3M[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,BdTR3C[&Created="Tue May 23 17:01:31 CEST 2017"]:1.5E-5)[&"Consensus support(%)"=66.0]:4.9E-5)[&"Consensus support(%)"=100.0]:0.029073,(Koz-2[&Created="Tue May 23 17:01:31 CEST 2017"]:4.0E-6,Koz-1[&Created="Tue May 23 17:01:31 CEST 2017"]:2.9E-5)[&"Consensus support(%)"=100.0]:0.018384)[&"Consensus support(%)"=100.0]:0.033958)[&"Consensus support(%)"=100.0]:0.01164)[&"Consensus support(%)"=73.0]:0.006661)[&"Consensus support(%)"=65.0]:0.008072)[&"Consensus support(%)"=99.0]:0.015779)[&"Consensus support(%)"=99.0]:0.013677)[&"Consensus support(%)"=98.0]:0.006214)[&"Consensus support(%)"=99.0]:0.031176)[&"Consensus support(%)"=100.0]:0.306644,(Bhyb36[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,BhybsixD[&Created="Tue May 23 17:01:31 CEST 2017"]:0.016639)[&!rotate=true,"Consensus support(%)"=100.0]:0.06122)[&!rotate=true,"Consensus support(%)"=66.0]:0.040434)[&!rotate=true,"Consensus support(%)"=100.0]:0.115179)[&!rotate=true,"Consensus support(%)"=100.0]:0.094708)[&!rotate=true,"Consensus support(%)"=99.0]:0.021095)[&!rotate=true,"Consensus support(%)"=94.0]:0.021012)[&!rotate=true,"Consensus support(%)"=100.0]:0.105443)[&!rotate=true,"Consensus support(%)"=72.0]:0.018588)[&"Consensus support(%)"=94.0]:0.04094,(ABR9_r[&Created="Tue May 23 17:01:31 CEST 2017"]:0.064347,((Geo-G30i2[&Created="Tue May 23 17:01:31 CEST 2017"]:0.011711,(Geo-G33i6[&Created="Tue May 23 17:01:31 CEST 2017"]:0.002199,Geo-G33i4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=100.0]:0.0577)[&"Consensus support(%)"=81.0]:0.057824,((Alb-AL2D[&Created="Tue May 23 17:01:31 CEST 2017"]:0.029142,('Alb-AL1A.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.048009,(Alb-AL2F[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0,Alb-AL2E[&Created="Tue May 23 17:01:31 CEST 2017"]:0.0)[&"Consensus support(%)"=100.0]:0.054256)[&"Consensus support(%)"=80.0]:0.01564)[&"Consensus support(%)"=100.0]:0.049368,('Ita-Sic-SLZ2.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.042576,(Ita-Sic-BNT3[&Created="Tue May 23 17:01:31 CEST 2017"]:0.042196,((Bd1-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.052086,(Tek-4[&Created="Tue May 23 17:01:31 CEST 2017"]:0.04957,(Tek-9[&Created="Tue May 23 17:01:31 CEST 2017"]:0.051,((BdTR7a[&Created="Tue May 23 17:01:31 CEST 2017"]:0.022507,(Tek-2[&Created="Tue May 23 17:01:31 CEST 2017"]:2.5E-5,(BdTR8i[&Created="Tue May 23 17:01:31 CEST 2017"]:0.001321,BdTRi[&Created="Tue May 23 17:01:31 CEST 2017"]:5.36E-4)[&"Consensus support(%)"=98.0]:0.007104)[&"Consensus support(%)"=100.0]:0.084624)[&"Consensus support(%)"=96.0]:0.041099,(Bd29-1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.041386,Ukr-Nvk1[&Created="Tue May 23 17:01:31 CEST 2017"]:0.029757)[&"Consensus support(%)"=48.0]:0.014976)[&"Consensus support(%)"=63.0]:0.017686)[&"Consensus support(%)"=39.0]:0.011026)[&"Consensus support(%)"=65.0]:0.014011)[&"Consensus support(%)"=98.0]:0.029201,(Ita-Sic-BNT8[&Created="Tue May 23 17:01:31 CEST 2017"]:0.045448,'Ita-Sic-BNT4.1'[&Created="Tue May 23 17:01:31 CEST 2017"]:0.017138)[&"Consensus support(%)"=69.0]:0.01896)[&"Consensus support(%)"=46.0]:0.011282)[&"Consensus support(%)"=71.0]:0.01955)[&"Consensus support(%)"=41.0]:0.017403)[&"Consensus support(%)"=47.0]:0.014715)[&"Consensus support(%)"=88.0]:0.02317)[&"Consensus support(%)"=55.00000000000001]:0.018931)[&"Consensus support(%)"=100.0]:0.109801);"""
regex = re.compile(".*?\[(.*?)\]")
results = re.findall(regex, tree_raw)
tree = tree_raw.replace('[','').replace(']','').replace("'",'')
if 0:
    for result in results:
        tree = tree.replace(result,'')
    print tree
    #print tree#tree_raw.replace("""[&Created="Tue May 23 17:01:31 CEST 2017"]""",'').replace("""[&"Consensus support(%)"=100.0]""",'')
    df = pd.read_csv('merged_table_w_depth_deDuplicated_removedDupCols.tr.csv')
    #print df['project_id']
    #print Tree(tree)
    #print traceback
    strain = list(df['strain'])
    protDict = []
    noInclude = []
    valL = []
    for i,val in enumerate(df['joshua_assembly_name']):
        if str(val) != 'nan':
            #print val
            try:
                ID = val.split('_')[1]
                IDnew = traceback[ID]
                if IDnew in list110:
                    #print i,val
                    df.set_value(i,'TreeName','Bdist%s'%(IDnew))
                    df.set_value(i,'Fake ID',IDnew)
                    tree = tree.replace(strain[i],'Bdist%s'%(IDnew))
                    noInclude.append(i)
                    protDict.append('%s     %s'%(IDnew,'Bdist'+IDnew))
            except:
                valL.append((i,val))
    #print len(df[df['TreeName'] != 'a'])
    sL = len(df[df['TreeName'] != 'a'])
    for i,val in enumerate(strain):
        if i not in noInclude:
            try:
                IDnew = oldDist[val.lower()]
                print IDnew + 'num1'
                if IDnew in list110:
                    df.set_value(i,'TreeName','Bdist%s'%(IDnew))
                    df.set_value(i,'Fake ID',IDnew)
                    tree = tree.replace(strain[i], 'Bdist%s' % (IDnew))
                    noInclude.append(i)
                else:
                    print IDnew + ' in bad'
            except:
                try:
                    IDnew = more[val.lower()]
                    print IDnew + 'num2'
                    if IDnew in list110:
                        df.set_value(i, 'TreeName', 'Bdist%s' % (IDnew))
                        df.set_value(i, 'Fake ID', IDnew)
                        tree = tree.replace(strain[i], 'Bdist%s' % (IDnew))
                        noInclude.append(i)
                except:
                    print val
                    valL.append((i, val))
    treenames  =list(df['TreeName'])
    projId = list(df['project_id'])
    """
    for i,val in enumerate(projId):
        if i not in noInclude:
            try:
                IDnew = oldDist[val.lower()]
                df.set_value(i,'TreeName','Bdist%s'%(IDnew))
                df.set_value(i,'Fake ID',IDnew)
            except:
                print val
                """
    #print Tree(tree)
    print more
    print traceback
    print oldDist
    print len(df[df['TreeName'] != 'a']) -  sL
    print len(valL)
    print len(list110)
    #print len(df[df['TreeName'] != 'a'])
    #print len(df[df['species']=='distachyon'])
    print len(set(bad))
    print len(df[df['TreeName'] == 'a'][df['species']=='distachyon'])
    #print df[df['species'] == 'distachyon'][df['TreeName'] == 'a']
    #print Tree(tree)
    print np.array(valL)
    #tree = tree.replace('Bd3-1_r','Bdistachyon')
    tree = tree.replace('Bdistachyon','Bdist314')
    t = Tree(tree)
    #print t
    treenames = [treename for treename in list(df['TreeName']) if treename != 'a']
    t.prune([node.name for node in t.traverse('postorder') if node.name in treenames and node.name] + ['Bdist314'])
    #t.prune([node.name for node in list(t.traverse('postorder')) if node.name])
    print t
    #print [node.name for node in list(t.traverse('postorder')) if node.name == 'Bdistachyon']
    #print len(list(t.traverse('postorder')))
    df.to_csv('correspondence.csv')
    with open('prot_dict','w') as f:
        f.write('\n'.join(protDict+['314     Bdist314']))
    #t.prune([node.name for node in t.traverse('postorder') if node.name in [] and node.name])
    t.write(format=1,outfile='grassnewick.nh')
    print bad
    print list110
    # make prot_dict
tree = Tree(open('grassnewick_091117temp7.nh','r').read())
print tree