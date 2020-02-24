library(data.table)
library(dplyr)

setwd("~/data/DroAr_data/")
#file with RPKM normalized read counts table
read_counts_normalized <- "Results/read_counts_RPKM_normalized.tsv"

#load read counts table
read_counts_normalized.dt <- fread(read_counts_normalized)

#classification of genes by time classes
#select ORFs
read_counts_normalized.dt_orfs <- read_counts_normalized.dt %>% 
  filter(!grepl("igr", feature))
#calculate logFC values between neighbor time points
logfc_df <- data.table(feature = counts.rpkms$feature,
                       logfc_40vs90 = log10(counts.rpkms$`90`)-log10(counts.rpkms$`40`),
                       logfc_90vs140 = log10(counts.rpkms$`140`)-log10(counts.rpkms$`90`),
                       logfc_140vs190 = log10(counts.rpkms$`190`)-log10(counts.rpkms$`140`))

logFC.dt <- data.table(feature = read_counts_normalized.dt_orfs$feature,
                       logfc_40vs90 = with(read_counts_normalized.dt_orfs, log10(t90)-log10(t40)),
                       logfc_90vs140 = with(read_counts_normalized.dt_orfs, log10(t140)-log10(t90)),
                       logfc_90vs190 = with(read_counts_normalized.dt_orfs, log10(t190)-log10(t90)),
                       logfc_140vs190 = with(read_counts_normalized.dt_orfs, log10(t190)-log10(t140)))

#classify ORFs according to the differences between different times of infection
putative_early_genes <- logFC.dt[logfc_90vs140 < 0 & logfc_140vs190]$feature
putative_middle_genes <- logFC.dt[logfc_90vs190 > 0 & logfc_90vs190 <= 1]$feature
putative_late_genes <- logFC.dt[logfc_90vs190 > 1]$feature

read_counts_normalized.dt_orfs.time_classes <- read_counts_normalized.dt_orfs %>% 
  mutate(time_class = ifelse(feature %in% putative_early_genes, "early", 
                             ifelse(feature %in% putative_middle_genes, "middle",
                                    ifelse(feature %in% putative_late_genes, "late", "ND"))))

fwrite(read_counts_normalized.dt_orfs.time_classes, "Results/gene_classes_and_RPKMs.tsv", sep = "\t", quote = F, row.names = F)

