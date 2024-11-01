# faecalmiRseq
Pipeline for identifying miRNAs in small RNA sequencing data of mouse faecal samples

This pipeline requires the following packages:  
cutadapt/1.8  
bowtie/1.1.1  
samtools/1.4  
subread/1.6.0  

And the following input files:  
* demultiplexed small RNA sequencing data in compressed fastq format (suffix fastq.gz)
* tRNA reference file e.g. from GtRNAdb 2.0 --> http://gtrnadb.ucsc.edu  
* rRNA reference e.g. from RNA central --> https://rnacentral.org  
* mouse genome reference e.g. from --> https://www.ncbi.nlm.nih.gov  
* mouse miRNA genome annotation file (mmu.gff3) e.g. from --> https://www.mirbase.org/download/  

all reference files (not including the miRNA genome annotation file) need to be indexed prior to use using bowtie-build  
e.g. bowtie-build <path_to_ref_fna_file> <path_to_output_file>  
output file should have the suffix .ind 

Please note that there may be newer versions of these packages that may contain bug fixes.

An R script for performing differential expression analysis of faecal miRNAs across conditions is also provided.
This R script requires the prior install of DESEq2, dplyr and tidyverse to the R environment.


