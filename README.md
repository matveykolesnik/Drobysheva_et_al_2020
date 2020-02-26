# phi14-2 transcriptome
Analysis of transcriptome of _Cellulophaga baltica_ cells infected with phi14:2 phage

_C. baltica_ phage phi14:2 belongs to a recently discovered family of crAss-like phages. phi14:2 phage possesses its own virion-packaged single-subunit RNA polymerase.
This repository contains a set of bash and R scripts that were used for the analysis of data obtained from high-throughput RNA sequencing of total RNA samples extracted from _C. baltica_ cells infected with phi14:2 phage at different times of post-infection.
These scripts estimate abundances of phi14:2 gene transcripts, classify these genes by the levels of expression at different stages of infection and build corresponding plots and heatmaps.

## Process_raw_data.sh
Bash script that performs QC of RNA-Seq reads before and after trimming and removes adapters and low-quality reads segments. Then script maps trimmed read onto reference sequences, prepares sorted and indexed BAM files and counts the numbers of mapped reads in each library.
**Requirements** fastqc, trimmomatic, bowtie2, samtools
**Input:** Raw reads files (FASTQ), reference sequences (FASTA)
**Output:** FastQC reports,sorted and indexed BAM files, table with the numbers of mapped reads per sample

## build_read_counts_tables.R
R script that counts reads overlapping with annotated phi14:2 genome features and performs RPKM normalization.
**Requirements:** R with Rsubread, edgeR, data.table and reshape2 packages
**Input:** Sorted and indexed BAM files, phi14:2 genomic features coordinates in GFF3 format, table with the numbers of mapped reads per sample
**Output:** Tables with reads counts, table with RPKM values

## ORFs_time_classes.R
R script that classifies phi14:2 genome features according to the transcript abundances at different time points
**Requirements:** R with data.table and dplyr packages
**Input:** Table with RPKM values
**Output:** Table with ORFs with the assigned stages of expression ("early", "middle" or "late")

## Draw_transcripts_abundances_plots.R
R script that draws plots and heatmaps
**Requirements:** R with data.table, dplyr, reshape2, ggplot2, Gviz and rtracklayer packages
**Input:** phi14:2 genomic features coordinates in GFF3 format, tables with RPKM values and ORFs time classes
**Output:** Plots and heatmaps
