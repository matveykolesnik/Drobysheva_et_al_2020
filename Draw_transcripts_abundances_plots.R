library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)

setwd("~/data/Drobysheva_et_al_2020/")
#File with normalized transcript abundances and time classes
classified_ORFs <- "Results/gene_classes_and_RPKMs.tsv"
#Annotation file
phi14_2_annotation <- "annotation/NC_021806_annotation_by_EVK.gff3"
#transcript abundances table
tr_ab_table <- "Results/read_counts_RPKM_normalized.tsv"

#load data
classified_ORFs_Rif_minus.dt <- fread(classified_ORFs) %>% 
  select(-c("without_phage", "t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif")) #remove uninfected sample and Rif+ samples
#g107 ORF is located within middle genes segment and we suspect that this gene is also should belong to middle gene class
classified_ORFs_Rif_minus.dt[feature == "g107"]$time_class <- "middle"
#normalize transcript abundances for each ORF to max. transcript abundance across all time points
classified_ORFs_Rif_minus.dt_maxnorm <- classified_ORFs_Rif_minus.dt
classified_ORFs_Rif_minus.dt_maxnorm[,c("t40", "t90", "t140", "t190")] <- classified_ORFs_Rif_minus.dt_maxnorm[,c("t40", "t90", "t140", "t190")]/apply(classified_ORFs_Rif_minus.dt_maxnorm[,c("t40", "t90", "t140", "t190")], 1, max)

classified_ORFs_Rif_minus.dt_maxnorm.melted <- melt(classified_ORFs_Rif_minus.dt_maxnorm, id.vars = c("feature", "time_class"), variable.name = "time_point", value.name = "abundance") %>% 
  mutate(time_class = factor(time_class, levels = c("early", "middle", "late")))
#Draw plot
tr_abundances <- ggplot(classified_ORFs_Rif_minus.dt_maxnorm.melted, aes(x = time_point, y = abundance, group = feature, colour = feature)) +
  geom_line() +
  facet_grid(cols = vars(time_class)) +
  theme(legend.position = "none") +
  xlab(label = "Time post-infection, min") +
  ylab(label = "Transcript abundances, % of max. expression") +
  scale_x_discrete(labels=c("40", "90", "140", "190")) +
  theme_bw() +
  theme(legend.position="none")

ggsave("Results/Pictures/Transcript_abundances_Rif_minus.svg", tr_abundances)

#Draw plot of transcript abundances in Rif+ samples normalized to local max. values in Rif- samples
classified_ORFs.dt <- fread(classified_ORFs) %>% 
  select(-without_phage)
classified_ORFs.dt[feature == "g107"]$time_class <- "middle"

classified_ORFs.dt_Rif_plus_norm_to_Rif_minus <- classified_ORFs.dt
classified_ORFs.dt_Rif_plus_norm_to_Rif_minus[,c("t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif")] <- 
  classified_ORFs.dt_Rif_plus_norm_to_Rif_minus[,c("t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif")]/apply(classified_ORFs.dt_Rif_plus_norm_to_Rif_minus[,c("t40", "t90", "t140", "t190")], 1, max)
classified_ORFs.dt_Rif_plus_norm_to_Rif_minus <- classified_ORFs.dt_Rif_plus_norm_to_Rif_minus %>% 
  select(-c("t40", "t90", "t140", "t190"))

classified_ORFs.dt_Rif_plus_norm_to_Rif_minus.melted <- melt(classified_ORFs.dt_Rif_plus_norm_to_Rif_minus, id.vars = c("feature", "time_class"), variable.name = "time_point", value.name = "abundance")%>% 
  mutate(time_class = factor(time_class, levels = c("early", "middle", "late")))

tr_abundances_rifplus <- ggplot(classified_ORFs.dt_Rif_plus_norm_to_Rif_minus.melted, aes(x = time_point, y = abundance, group = feature)) +
  geom_line(aes(colour=time_class, alpha=feature), show.legend = FALSE) +
  scale_colour_manual(values = c(early="darkgreen", middle="purple", late="blue", guide="none")) +
  facet_wrap(~time_class, ncol=3) +
  scale_x_discrete(labels=c("40", "90", "140", "190")) +
  xlab(label = "Time after infection, min") +
  ylab(label = "Transcript abundance, %")

ggsave("Results/Pictures/Transcript_abundances_Rif_plus_norm_to_Rif_minus.svg", tr_abundances_rifplus)

library(Gviz)
library(rtracklayer)
options(ucscChromosomeNames = F)

tr_ab_table.dt <- fread(tr_ab_table) %>% 
  select(-c("without_phage", "t40_Rif", "t90_Rif", "t140_Rif", "t190_Rif")) #remove uninfected sample and Rif+ samples

tr_ab_table.dt_localmaxnorm <- tr_ab_table.dt
tr_ab_table.dt_localmaxnorm[,c("t40", "t90", "t140", "t190")] <- tr_ab_table.dt_localmaxnorm[,c("t40", "t90", "t140", "t190")]/apply(tr_ab_table.dt_localmaxnorm[,c("t40", "t90", "t140", "t190")], 1, max)
#Draw heatmaps

phi14_2_annotation.GR <- readGFFAsGRanges(phi14_2_annotation)
phi14_2_annotation.track <- AnnotationTrack(start = phi14_2_annotation.GR@ranges@start, 
                                            width = phi14_2_annotation.GR@ranges@width,
                                            chromosome = "NC_021806.1",
                                            strand = as.character(phi14_2_annotation.GR@strand),
                                            id = phi14_2_annotation.GR$ID,
                                            stacking="dense",
                                            fill = "grey",
                                            col = "black",
                                            lex = 10,
                                            lty = 1,
                                            shape = "fixedArrow",
                                            arrowHeadWidth = 100, 
                                            lwd = 1,
                                            fontcolor.feature="black",
                                            rotation.item = 90,
                                            cex = 2,
                                            labelPos = "below")

axisTrack <- GenomeAxisTrack(cex = 2, labelPos="below", size = 1)
hm_data <- GRanges(seqnames = phi14_2_annotation.GR@seqnames, ranges = phi14_2_annotation.GR@ranges, strand = rep("*", length(phi14_2_annotation.GR@ranges)))
hm_data$feature <- phi14_2_annotation.GR$ID
#Draw heatmap normalized to local max.
values(hm_data) <- merge(hm_data, tr_ab_table.dt_localmaxnorm, by.x = "feature", sort = F)[c("t40", "t90", "t140", "t190")]

dTrack <- DataTrack(hm_data, 
                    name = "rpkm/local_max(rpkms)",
                    type = "heatmap",
                    showSampleNames = T)

svg("Results/Pictures/transcript_abundances_heatmap_norm_on_norm_on_local_max.svg", width = 150, height = 40)
plotTracks(c(phi14_2_annotation.track, axisTrack, dTrack), featureAnnotation = "id")
dev.off()

#maximum value across transcript abundances of all ORFs across all time points
global_max_value <- max(apply(tr_ab_table.dt[,c("t40", "t90", "t140", "t190")], 1, max))
tr_ab_table.dt_globalmaxnorm <- tr_ab_table.dt
tr_ab_table.dt_globalmaxnorm[,c("t40", "t90", "t140", "t190")] <- tr_ab_table.dt_globalmaxnorm[,c("t40", "t90", "t140", "t190")]/global_max_value

hm_data <- GRanges(seqnames = phi14_2_annotation.GR@seqnames, ranges = phi14_2_annotation.GR@ranges, strand = rep("*", length(phi14_2_annotation.GR@ranges)))
hm_data$feature <- phi14_2_annotation.GR$ID
#Draw heatmap normalized to global max.
values(hm_data) <- merge(hm_data, tr_ab_table.dt_globalmaxnorm, by.x = "feature", sort = F)[c("t40", "t90", "t140", "t190")]

dTrack <- DataTrack(hm_data, 
                    name = "rpkm/local_max(rpkms)",
                    type = "heatmap",
                    showSampleNames = T)

svg("Results/Pictures/transcript_abundances_heatmap_norm_on_norm_on_global_max.svg", width = 150, height = 40)
plotTracks(c(phi14_2_annotation.track, axisTrack, dTrack), featureAnnotation = "id")
dev.off()
