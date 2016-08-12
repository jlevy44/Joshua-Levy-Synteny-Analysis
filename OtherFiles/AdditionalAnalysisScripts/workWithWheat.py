from pyfaidx import Fasta
"""
wheat = Fasta('Taestivum_296_v2.softmasked.fa')

a=1
for key in wheat.keys():
    if '7al' in key:
        print key, wheat[key][0:150]
"""

# lets try to convert wheat to 7 chromosomes
chromosomesStartEndCoords = {1:[],2:[],3:[],4:[],5:[],6:[],7:[]}
open('karyotype.Taestivum_296.txt','w').close()
karyOut = open('karyotype.Taestivum_296.txt','w')

for key in chromosomesStartEndCoords.keys():
    inputFai = open('Taestivum_296_v2.softmasked.fa.fai','r')
    readLast = 0
    for line in inputFai:
        lineList = line.split()
        if str(key)+'al' in lineList[0].split('_')[2] and chromosomesStartEndCoords[key] == []:
            chromosomesStartEndCoords[key].append(int(lineList[2]))
            chromosomesStartEndCoords[key].append(0)
        if readLast ==1 and str(key)+'ds' in lineList[0].split('_')[2]:
            if int(lineList[2])> chromosomesStartEndCoords[key][1]:
                chromosomesStartEndCoords[key][1]=int(lineList[2])
        elif readLast ==1:
            break
        if readLast == 0 and str(key)+'ds' in lineList[0].split('_')[2]:
            readLast = 1
            chromosomesStartEndCoords[key][1]=int(lineList[2])
    inputFai.seek(0)
    inputFai.close()
    # chr - Chr1 Chr1 0 30427671 chr1
    karyOut.write('chr - ta%d ta%d %d %d chr%d\n'%(key,key,chromosomesStartEndCoords[key][0],
                                                 chromosomesStartEndCoords[key][1],key))
karyOut.close()




inputFai = open('Taestivum_296_v2.softmasked.fa.fai','r')

bedFileInput = open('PAC2_0.283-PAC2_0.296_5.bed','r')
open('PAC2_0.283-PAC2_0.296_5.link.txt','w').close()
linkFileOut = open('PAC2_0.283-PAC2_0.296_5.link.txt','w')
for line in bedFileInput:
    lineList = line.split()
    lineList2 = lineList[-1].split('-')

    for lineFai in inputFai:
        if lineList2[1] in lineFai:
            outCoords = [int(lineFai.split()[2])+int(lineList2[2]),int(lineFai.split()[2])+int(lineList2[3])]
            inputFai.seek(0)
            break
    outputTuple = (lineList[0].split('-')[1], lineList[1], lineList[2], lineList2[1].split('_')[2][0], outCoords[0],
                   outCoords[1])
    linkFileOut.write('%s %s %s ta%s %d %d\n' % outputTuple)

linkFileOut.close()
inputFai.close()
bedFileInput.close()

print chromosomesStartEndCoords

"""
283-Bd1	35933291	35964024	296-ta_iwgsc_5as_v1_1552572-13295-44443
283-Bd1	16613913	16643218	296-ta_iwgsc_2ds_v1_5318666-1888-34988
283-Bd1	35940758	35957021	296-ta_iwgsc_7al_v1_4559370-2730-17574
283-Bd3	56364839	56382958	296-ta_iwgsc_6al_v1_5759130-3349-36390
283-Bd1	35952306	35964024	296-ta_iwgsc_1bs_v1_3483406-6717-18362
283-Bd5	20186823	20210770	296-ta_iwgsc_2al_v1_6438836-4411-34205
283-Bd1	35933291	35942931	296-ta_iwgsc_2as_v1_5199852-7048-17282
283-Bd1	35950411	35961554	296-ta_iwgsc_4bs_v1_4866414-3176-13526
283-Bd1	35936638	35944690	296-ta_iwgsc_5bl_v1_10926321-1223-9176
283-Bd1	66539276	66567922	296-ta_iwgsc_4as_v2_5993421-6036-36536
283-Bd1	25451452	25474694	296-ta_iwgsc_7dl_v1_3311257-448-16245
283-Bd1	35942241	35957021	296-ta_iwgsc_3dl_v1_6929075-2084-15979
283-Bd1	3910963	3928271	296-ta_iwgsc_5as_v1_1552572-22819-44443
283-Bd1	35933291	35942931	296-ta_iwgsc_1bs_v1_3430402-834-11068
283-Bd1	7898748	7949977	296-ta_iwgsc_4bs_v1_4875646-4648-24523
283-Bd1	11627943	11668367	296-ta_iwgsc_4ds_v1_2311530-3029-30515
283-Bd3	3273824	3287204	296-ta_iwgsc_6as_v1_4338458-2864-41960
283-Bd1	35957854	35964024	296-ta_iwgsc_5ds_v1_2780529-285-6484
283-Bd1	11111882	11126532	296-ta_iwgsc_5bl_v1_10732661-4219-20247
283-Bd1	69461389	69493462	296-ta_iwgsc_4dl_v3_14370381-8421-25912
283-Bd1	35940758	35954151	296-ta_iwgsc_6bs_v1_2926507-807-11317
283-Bd1	35952306	35961554	296-ta_iwgsc_1ds_v1_1902000-187-9384
283-Bd1	70402950	70446047	296-ta_iwgsc_4as_v2_5959767-538-32015
283-Bd1	35950411	35958874	296-ta_iwgsc_4bl_v1_7012057-179-8281
283-Bd2	11500054	11509154	296-ta_iwgsc_2as_v1_5199852-7048-16115
283-Bd5	22407366	22423788	296-ta_iwgsc_2bl_v1_8045172-9207-31179
283-Bd1	19458599	19475896	296-ta_iwgsc_2as_v1_5222718-4238-24001
283-Bd3	12787260	12806743	296-ta_iwgsc_7al_v1_4512428-1165-21657
283-Bd1	35938638	35944690	296-ta_iwgsc_7as_v1_4202288-497-6252
283-Bd1	4615125	4628565	296-ta_iwgsc_5dl_v1_4528289-3-15287
283-Bd1	7155099	7179219	296-ta_iwgsc_4bs_v1_4942770-3847-28275
283-Bd5	19813575	19828142	296-ta_iwgsc_2bl_v1_7949001-4751-17510
283-Bd5	26576733	26588410	296-ta_iwgsc_2dl_v1_9907310-2350-21850
283-Bd1	3905437	3908068	296-ta_iwgsc_6bs_v1_3019597-4033-6676
283-Bd4	35922940	35940436	296-ta_iwgsc_5bl_v1_10821942-2956-18936
283-Bd1	74737505	74749251	296-ta_iwgsc_4dl_v3_14467290-2-16868
283-Bd5	15332334	15366571	296-ta_iwgsc_2bl_v1_8090499-6600-27781
283-Bd1	9959592	9980794	296-ta_iwgsc_4ds_v1_2312826-2191-24461
283-Bd1	65080370	65102709	296-ta_iwgsc_4as_v2_5984737-7744-29137
283-Bd5	26886772	26902310	296-ta_iwgsc_2al_v1_6377650-1007-20654
283-Bd3	51009521	51026250	296-ta_iwgsc_6dl_v1_3290973-4819-21877
283-Bd4	36285254	36303204	296-ta_iwgsc_5bl_v1_10796095-1738-16940
283-Bd1	3899742	3911662	296-ta_iwgsc_6as_v1_4428294-6690-13087
283-Bd2	11497748	11505220	296-ta_iwgsc_5bl_v1_10926321-1223-9176
283-Bd3	56367381	56377308	296-ta_iwgsc_6dl_v1_3324842-0-11398
283-Bd1	15595115	15608645	296-ta_iwgsc_2bs_v1_5202128-4706-19320
283-Bd2	11497748	11509154	296-ta_iwgsc_5as_v1_1552572-13295-24706
283-Bd1	74845239	74861674	296-ta_iwgsc_4dl_v3_14404266-582-20304
283-Bd1	35950411	35958874	296-ta_iwgsc_4as_v2_5949307-3176-11170
283-Bd1	35940758	35947385	296-ta_iwgsc_4ds_v1_2300695-416-6942
283-Bd1	35940758	35950196	296-ta_iwgsc_7bs_v1_3168908-0-8782
283-Bd1	35950411	35957766	296-ta_iwgsc_6as_v1_4428294-6690-13087
283-Bd1	35958514	35964024	296-ta_iwgsc_2al_v1_6336037-0-5397
283-Bd2	11500054	11509154	296-ta_iwgsc_1bs_v1_3430402-2001-11068
283-Bd4	37728933	37745738	296-ta_iwgsc_5bl_v1_10794259-3927-20207
283-Bd5	24108915	24129711	296-ta_iwgsc_2bl_v1_7990526-3922-22959
283-Bd4	34099042	34131601	296-ta_iwgsc_5dl_v1_4496039-2246-28478
283-Bd1	35958514	35964024	296-ta_iwgsc_7al_v1_4535756-5823-10815
283-Bd1	3775447	3787205	296-ta_iwgsc_5bl_v1_10845579-60-12653
283-Bd1	6188077	6211105	296-ta_iwgsc_5dl_v1_4580238-7519-25467
"""