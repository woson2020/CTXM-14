#-*-coding:utf-8 -*-
#2018-11-19
#Pacbio_sequence analysis yang, find illumina_sequencing umi
from sys import argv
import os
import regex

input1,input2,output1=argv[1:]
out1=open(output1,"w")
out2=open(output2,"w")

# 1 base error is acceptable in barcode calling 
f_re="(?e)(?P<p1>TAGGACAAATCCGCCTTAATTAA){e<=1}(?P<umi_f>[ATCG]{20})(?P<p2>CTGCAGGCATGCATTTAAATCTC){e<=1}"
r_re="(?e)(?P<p3>GAGATTTAAATGCATGCCTGCAG){e<=1}(?P<umi_r>[ATCG]{20})(?P<p4>TTAATTAAGGCGGATTTGTCCTA){e<=1}"
f_re_alter="(?e)(?P<p5>TAGGACAAATCCGCC){e<=1}(?P<umi_f_alter>[ATCG]{20})(?P<p6>CTGCAGGCATGCCTCGAGATG){e<=1}"#for manually constructed genotypes
r_re_alter="(?e)(?P<p7>CATCTCGAGGCATGCCTGCAG){e<=1}(?P<umi_r_alter>[ATCG]{20})(?P<p8>GGCGGATTTGTCCTA){e<=1}"#for manually constructed genotypes

mud_dict={"A":"T","T":"A","C":"G","G":"C","N":"N"}
reverComple=lambda s: "".join([mud_dict[c] for c in s])[::-1]


with open(input1,"r") as r1,open(input2,"r") as r2:
    ilumi_umi_set=set()
    for i,j in enumerate(zip(r1,r2)):      
        if i % 4==1 :
            j0=j[0].strip()
            j1=j[1].strip()
            f_pre=regex.search(f_re,j0,flags=0)
            r_pre=regex.search(r_re,j1,flags=0)
            f_pre_alter=regex.search(f_re_alter,j0,flags=0)
            r_pre_alter=regex.search(r_re_alter,j1,flags=0)
            if f_pre is not None and r_pre is not None:
                f=f_pre.groupdict()["umi_f"]
                r=reverComple(r_pre.groupdict()["umi_r"])
                if f==r:
                    ilumi_umi_set.add(f)
                    out1.write(f+"\n")
            if f_pre is  None and r_pre is None and f_pre_alter is not None and r_pre_alter is not None:
                f_alter=f_pre_alter.groupdict()["umi_f_alter"]
                r_alter=reverComple(r_pre_alter.groupdict()["umi_r_alter"])
                if f_alter==r_alter:
                    out1.write(f_alter+"\n")
out1.close()


