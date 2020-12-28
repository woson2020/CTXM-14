#-*-coding:utf-8 -*-
from sys import argv
import os
import re
import time
import multiprocessing

"""
blasr can be download fron https://github.com/PacificBiosciences/blasr (suggest install via bioconda).
ccs can also be install via bioconda.
samtools with version=1.9 (using htslib 1.9) ,or may have some conflict
"""
raw_bam,ref_file,threads=argv[1:]#raw_bam is the raw file of PacBio-subreads(must end with .bam),this file should be in the same folder with this python script


os.system("blasr --bam --out %s.blasr.bam %s %s"%(raw_bam[:-4],raw_bam,ref_file))#mapping to the wild-type of gfp
os.system("samtools view -h -o %s.blasr.sam %s.blasr.bam"%(raw_bam[:-4],raw_bam[:-4]))#tranform bam file to sam file 
os.system("samtools view -h -o %s.sam %s"%(raw_bam[:-4],raw_bam))
p_strand=open("%s.p.sam"%(raw_bam[:-4]),"w")#positive strand
n_strand=open("%s.n.sam"%(raw_bam[:-4]),"w")#negative strand

with open("%s.sam"%(raw_bam[:-4]),"r") as sam_file:
    subreads_s_dict={}
    for read_line in sam_file:
        read_line=read_line.strip()
        if read_line[0]=="m":
            read_line1_list=read_line.split()
            subreads_s_dict[read_line1_list[0]]=read_line#build a dict (key is the ZMW name[the first column],value is subread)
        else:
            #output the header 
            p_strand.write(read_line+"\n")
            n_strand.write(read_line+"\n")

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])

with open("%s.blasr.sam"%(raw_bam[:-4]),"r") as blasr_sam:#把每一个zmw分成正负链两部分

        for sam_line in blasr_sam:
            if sam_line[0]=="m":
                sam_line_list=sam_line.strip().split()
                if sam_line_list[1]=="16":#positive strand
                    n_strand.write(subreads_s_dict[sam_line_list[0]]+"\n")   
                elif sam_line_list[1]=="0":#negative strand
                    p_strand.write(subreads_s_dict[sam_line_list[0]]+"\n")
p_strand.close()                 
n_strand.close()
#remove intermediate file
os.system("echo "" > *%s.blasr.sam|rm *%s.blasr.sam"%(raw_bam[:-4],raw_bam[:-4]))
os.system("echo "" > *%s.sam|rm *%s.sam"%(raw_bam[:-4],raw_bam[:-4]))

def run_ccs(name):
    os.system("samtools view -bS %s.sam -o %s.bam"%(name[:-4],name[:-4]))
    #plz pay attention to the number of thread
    os.system("ccs --min-length 900 --max-length 1600 -j %s --min-passes 5  %s.bam %s_ccs.bam"%(threads,name[:-4],name[:-4]))
    os.system("samtools view -h -o %s_ccs.sam %s_ccs.bam"%(name[:-4],name[:-4]))

#run ccs of positive strand set and negative strand set
p_n=[raw_bam[:-4]+".p.sam",raw_bam[:-4]+".n.sam"]
multiple_p=multiprocessing.Pool(2)
for file in p_n:
    ress=multiple_p.apply_async(run_ccs,args=(file,))
     
print('Waiting for all subprocesses done...')
multiple_p.close()
multiple_p.join()
print('All subprocesses done.')

#remove intermediate file
os.system("echo "" > *%s.n.bam|rm *%s.n.bam"%(raw_bam[:-4]),raw_bam[:-4])
os.system("echo "" > *%s.p.bam|rm *%s.p.bam"%(raw_bam[:-4],raw_bam[:-4]))
os.system("echo "" > *%s.n.sam|rm *%s.n.sam"%(raw_bam[:-4],raw_bam[:-4]))
os.system("echo "" > *%s.p.sam|rm *%s.p.sam"%(raw_bam[:-4],raw_bam[:-4]))

geno_pattern=re.compile("TGGGCAGCGCGCCGCTTTAT\w{795}CCAGGCATCAAATAAAACGA")
umi_pattern=re.compile("GACAAATCCGCCTTAATTAA\w{20}CTGCAGGCATGCATTTAAAT")
outfile=open("%s.final_result"%(raw_bam[:-4]),"w")
with open("%s.n_ccs.sam"%(raw_bam[:-4]),"r") as inputr:#nagetive strand ccs result
    for nega in inputr:    
        if nega[0]=="m":
            nega_line_list=nega.strip().split()
            nega_seq=nega_line_list[9]#the result sequence after ccs 
            zmw_nega=nega_line_list[0].split("/")[1]#the number of ZMW
            nega_reverse=Reverse_complemrnt(nega_seq)
            geno2=re.search(geno_pattern,nega_reverse,flags=0)
            umi2=re.search(umi_pattern,nega_reverse,flags=0)
            if geno2 is not None and umi2 is not None:
                geno3 = geno2.group()[20:815]
                umi3 = umi2.group()[20:40]
                outfile.write(">"+str(zmw_nega)+"n"+"\n"+"n_"+geno3+" "+umi3+"\n")
with open("%s.p_ccs.sam"%(raw_bam[:-4]),"r") as inputf:#positive strand ccs result      
    for posi in inputf:
        if posi[0]=="m":
            posi_line_list=posi.strip().split()
            posi_seq=posi_line_list[9]
            zmw_posi=posi_line_list[0].split("/")[1]
            geno=re.search(geno_pattern,posi_seq,flags=0)
            umi=re.search(umi_pattern,posi_seq,flags=0)
            if geno is not None and umi is not None:
                geno1 = geno.group()[20:815]
                umi1 = umi.group()[20:40]
                outfile.write(">"+str(zmw_posi)+"p"+"\n"+"p_"+geno1+" "+umi1+"\n")
                

os.system("rm *%s.*_ccs.bam.pbi"%(raw_bam[:-4]))
os.system("rm *%s.*_ccs.sam"%(raw_bam[:-4]))
outfile.close()