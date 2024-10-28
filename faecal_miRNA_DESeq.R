#####################################
## DESEQ analysis of faecal miRNAs ##
#####################################

# ensure that the columns of the counts file are in the same order as the condition variable
# i.e. the 5 samples in condition_A followed by B then C
# the counts file contains the raw counts output by featurecounts
library(DESeq2)
library(dplyr)
condition <- as.factor(c(rep("condition_A", 5), rep("condition_B", 5), rep("condition_C", 5)))
des<-formula(~ condition)
countsTable <-read.table('path_to_counts_file.txt', check.names=FALSE,row.names=1, header = TRUE)

# first we remove the duplicate loci
# ensure the column containing the MIMAT IDs is named microRNA_ID
# convert the row names to be the first column
countsTable <- data.frame(microRNA_ID  = row.names(countsTable), countsTable)
# search the MIMAT IDs and if it ends in _* then delete this
countsTable$microRNA_ID <- ifelse(grepl("MIMAT", countsTable$microRNA_ID), gsub("_.*", "", countsTable$microRNA_ID), countsTable$microRNA_ID)
# now are left with duplicate read IDs
# sum reads across all samples and then choose the duplicate loci with highest counts
# IMPORTANT - replace n1:n2 with column numbers to reflect the range of columns in the countsTable created above that contain the raw count data
countsTable$summedreads <- rowSums(countsTable[,n1:n2])
microRNAs_filt <- countsTable %>% group_by(microRNA_ID) %>% top_n(1, summedreads)
# if accessions have the same (highest number) of reads then remove one
microRNAs_filt_unique <- microRNAs_filt[!duplicated(microRNAs_filt$microRNA_ID), ]
# remove the excess columns (replace n1:n2 with relevant columns)
countsTable <- microRNAs_filt_unique[,n1:n2]
library(tidyverse)
countsTable <- countsTable %>% remove_rownames %>% column_to_rownames(var="microRNA_ID")

# perform DESEq, update n1:n2 to include the columns containing raw counts
myNames<-colnames(countsTable[c(n1:n2)])
colDataNames<-data.frame(row.names=myNames, condition=condition)
ddsHTSeq<-DESeqDataSetFromMatrix(countData = countsTable[c(n1:n2)],
                                 colData=colDataNames, design=des, ignoreRank = FALSE)
colDataNames
dds<-DESeq(ddsHTSeq,betaPrior=FALSE)
resultsNames(dds)

# write normalised counts
normCounts<-as.data.frame(counts(dds,normalized=TRUE))
# optionally save normCounts locally with: write.csv(normCounts, "deseq_normalised_counts.csv")

# make each comparison
res <- results(dds, contrast = c("condition", "condition_A", "condition_B"), cooksCutoff = TRUE, independentFiltering = TRUE)
summary(res)
head(res)
# repeat as desired
# optionally save locally with: write.csv(res, "deseq_results_condition_A_vs_condition_B.csv")