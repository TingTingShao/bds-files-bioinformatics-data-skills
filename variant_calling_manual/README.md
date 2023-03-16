# NGS variant calling workflow           
## Case          
Two exome sequencing samples of a T-Cell Acute Lymphoblastic Leukemia patient. One sample `TLE66_T` is of tumor, one healthy control sample `TLE66_N`
## Goal: process the raw fastq, call & annotate SNPs.
We want to ultimately find interesting variants that may be causal to the T-ALL.
# preparation
    > export PATH=/staging/leuven/stg_00079/teaching/miniconda3/envs/large_omics_2023_b/bin:$PATH
    # create a work folder!
    > cd $VSC_SCRATCH
    > mkdir -p variant_calling_manual
    > cd variant_calling_manual
    # We will be using the genome sequence & index a few times - 
    # So I'll use an environment variable
    DB=/staging/leuven/stg_00079/teaching/hg19_9/chr9.fa
    
