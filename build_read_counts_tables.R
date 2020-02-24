library(Rsubread)
library(data.table)
library(dplyr)
library(reshape2)
library(edgeR)
library(openxlsx)

#This script counts RNA-Seq reads mapped to the features of phi14:2 genome and performs RPKM normalization

setwd("~/data/DroAr_data/")

#directory with BAM files
aln_path <- "alignments/"
#file with number of mapped reads per sample
reads_per_lib <- "reads_per_lib.tsv"
#list of BAM files
bams <- list.files(aln_path, pattern = "\\.bam$")
#Annotation of phi14:2 phage in GFF3 format
phi14_2_annotation <- "annotation/NC_021806_annotation_by_EVK.gff3"

lib_sizes <- fread(reads_per_lib, sep="\t", header = F, col.names = c("libname", "size"))
lib_sizes <- setNames(data.table(t(lib_sizes[,-"libname"])), lib_sizes[["libname"]])

#Counting reads crossing ORFs and intergenic regions in phi14:2 genome
counts <- featureCounts(files = paste0(aln_path, bams),
                        annot.ext = phi14_2_annotation,
                        isGTFAnnotationFile = T,
                        GTF.featureType = "gene",
                        GTF.attrType = "ID",
                        allowMultiOverlap = T, #since ORFs are suspected to be organized in operons
                        strandSpecific = 2, #reads are reverse-complement
                        nthreads = 4)
counts.dge <- DGEList(counts = counts$counts,
                      genes = counts$annotation[, c("GeneID", "Length")])

#Save raw read counts table 
raw_counts_table <- as.data.table(counts.dge[["counts"]], keep.rownames = "feature")
colnames(raw_counts_table) <- c("feature", "without_phage", "t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif", "t40", "t90", "t140", "t190")
#write tables to file
write.table(x = raw_counts_table, file = "Results/raw_read_counts.tsv", sep = "\t", quote = F, row.names = F)
write.xlsx(x = raw_counts_table, file = "Results/raw_read_counts.xls")

#RPKM normalization
counts.dge[["counts"]][counts.dge[["counts"]] == 0] <- 1 #pseudocounts
#Normalize to the total number of mapped reads
counts.dge$samples$lib.size <- as.numeric(lib_sizes)
counts.rpkms <- rpkm(counts.dge,
                     gene.length = counts.dge$genes$Length)
counts.rpkms <- as.data.table(counts.rpkms, keep.rownames = "feature")
colnames(counts.rpkms) <- c("feature", "without_phage", "t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif", "t40", "t90", "t140", "t190")
#write tables to file
write.table(x = counts.rpkms, file = "Results/read_counts_RPKM_normalized.tsv", sep = "\t", quote = F, row.names = F)
write.xlsx(x = counts.rpkms, file = "Results/read_counts_RPKM_normalized.xls")
