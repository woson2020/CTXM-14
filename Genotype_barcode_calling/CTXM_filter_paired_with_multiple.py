#-*-coding:utf-8 -*-
from sys import argv
import numpy as np

ccs_result,output1=argv[1:]

wild_type="GCGCAGACGAGTGCGGTGCAGCAAAAGCTGGCGGCGCTGGAGAAAAGCAGCGGAGGGCGGCTGGGCGTCGCGCTCATCGATACCGCAGATAATACGCAGGTGCTTTATCGCGGTGATGAACGCTTTCCAATGTGCAGTACCAGTAAAGTTATGGCGGCCGCGGCGGTGCTTAAGCAGAGTGAAACGCAAAAGCAGCTGCTTAATCAGCCTGTCGAGATCAAGCCTGCCGATCTGGTTAACTACAATCCGATTGCCGAAAAACACGTCAACGGCACAATGACGCTGGCAGAACTGAGCGCGGCCGCGTTGCAGTACAGCGACAATACCGCCATGAACAAATTGATTGCCCAGCTCGGTGGCCCGGGAGGCGTGACGGCTTTTGCCCGCGCGATCGGCGATGAGACGTTTCGTCTGGATCGCACTGAACCTACGCTGAATACCGCCATTCCCGGCGACCCGAGAGACACCACCACGCCGCGGGCGATGGCGCAGACGTTGCGTCAGCTTACGCTGGGTCATGCGCTGGGCGAAACCCAGCGGGCGCAGTTGGTGACGTGGCTCAAAGGCAATACGACCGGCGCAGCCAGCATTCGGGCCGGCTTACCGACGTCGTGGACTGTGGGTGATAAGACCGGCAGCGGCGACTACGGCACCACCAATGATATTGCGGTGATCTGGCCGCAGGGTCGTGCGCCGCTGGTTCTGGTGACCTATTTTACCCAGCCGCAACAGAACGCAGAGAGCCGCCGCGATGTGCTGGCTTCAGCGGCGAGAATCATCGCCGAAGGGCTGTAA"
with open(ccs_result,"r") as ccs_file:#
    umi_list=[]
    ctx_list=set()
    umi_genot_dict={}#key is barcode,value is the list of genotype 
    multig=0
    for line in ccs_file:
        line=line.strip()
        if len(line)>700 :
            umi_g_list=line.split(" ")
            if len(umi_g_list[0])==795 and len(umi_g_list[1])==20:
                umi_pre=umi_g_list[1][:20]
                ctx_pre=umi_g_list[0][2::]
                umi_list.append(umi_pre)
                if umi_pre in umi_genot_dict.keys():#build barcode-genotype dict
                    umi_genot_dict[umi_pre].append(ctx_pre)
                    multig+=1
                else:
                    umi_genot_dict[umi_pre]=[ctx_pre] 
            else:
                print("error")

umi_mutant_dict={}#key is barcode,value is the unique genotype 

for i in umi_genot_dict.keys():
    if len(np.unique(umi_genot_dict[i]))==1:#barcode pair with unique genotyoe
        mutant=0
        call=[]
        genotype=""
        for g in range(795):
            if umi_genot_dict[i][0][g]!=wild_type[g]:
                s=wild_type[g]+str(g+1)+umi_genot_dict[i][0][g]
                genotype=genotype+s+" "
                mutant+=1
        if genotype=="":
            genotype="WT"

        genotype=">"+str(mutant)+"> "+genotype
        # if mutant>=1:
        umi_mutant_dict[i]=genotype
    else:#barcode pair with multiple genotyoe
        geno_dict={}
        uni_gen=np.unique(umi_genot_dict[i])
        for ii in uni_gen:
            geno_dict[ii]=umi_genot_dict[i].count(ii)
        sort_geno_dict=sorted(geno_dict.items(),key=lambda geno_dict:geno_dict[1],reverse=True)
        major_genotype=sort_geno_dict[0][1]#genotype with most reads
        if major_genotype>sort_geno_dict[1][1]:#major_genotype have the maximum reads
            mutant1=0
            genotype1=""
            for g1 in range(795):
                if sort_geno_dict[0][0][g1]!=wild_type[g1]:# call mutants
                    s1=wild_type[g1]+str(g1+1)+sort_geno_dict[0][0][g1]
                    genotype1=genotype1+s1+" "
                    mutant1+=1
            if genotype1=="":
                genotype1="WT"
            genotype1=">"+str(mutant1)+"> "+genotype1
            umi_mutant_dict[i]=genotype1     

g_umi_dict={}#key is genotype,value is the list barcode

for um in umi_mutant_dict.keys():
    if umi_mutant_dict[um] in g_umi_dict.keys():
        g_umi_dict[umi_mutant_dict[um]].append(um)
    else:
        g_umi_dict[umi_mutant_dict[um]]=[um]
with open(output1,"w") as out2:
    for ggg in g_umi_dict.keys():
        out2.write(ggg+"*")
        for uuu in g_umi_dict[ggg]:
            out2.write(" "+uuu)
        out2.write("\n")

