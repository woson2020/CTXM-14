1) CTXM_PacBio_genotype_barcode_calling.py is used to call barcode and its corresponding genotype according to PacBio sequencing data. 
Since we have sequenced three PacBio SMRT cells of our mutant library, each cells must run with this script and then combine 
those three cells' result to generate a summary file.

2) CTXM_filter_paired_with_multiple.py is used to filter genotype when a barcode corresponds to multiple genotypes.
The summary file from 1) is the input of this script, the output is  genotype and all its unique barcode.

3)CTXM_illumina_barcode_calling.py is used to call barcode according to Illumina sequencing data (from competitive culture).
Through this script, the frequency of various barcode was got and then use to calculate realtive growth.

#######
Ref_CTXM.fasta is the wildtype sequence of CTXM-14, which is necessary for CTXM_PacBio_genotype_barcode_calling.py.