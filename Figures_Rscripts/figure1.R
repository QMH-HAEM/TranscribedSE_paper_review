# R codes for generating Figure 1 in the transcribed SE paper
library(edgeR)
library(ggplot2)
library(ggsignif)
library(hciR)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
options(scipen = 9)
options(timeout=6000)

# Define general theme for plotting
theme01 <- theme(axis.text = element_text(size = 22),
                 axis.title = element_text(size = 24),
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 14),
                 legend.position = "none")

# Comparison between HK-AML, BeatAML and TCGA-LAML
# Load raw SE count matrices for the 3 cohorts

# Perform normalization
dge.hkaml <- DGEList(hkaml_count, genes = rownames(hkaml_count), lib.size = hkaml_lib.size)
dge.beataml <- DGEList(beataml_count, genes = rownames(beataml_count), lib.size = beataml_lib.size$lib.size)
dge.tcga <- DGEList(tcga_count, genes = rownames(tcga_count), lib.size = tcga_lib.size$lib.size)

dge.hkaml <- cpm(dge.hkaml)
dge.beataml <- cpm(dge.beataml)
dge.tcga <- cpm(dge.tcga)

hkaml_sum <- colSums(dge.hkaml)
beataml_sum <- colSums(dge.beataml)
tcga_sum <- colSums(dge.tcga)

# Compile dataframe for plotting
df_tmm <- data.frame(dataset = c(rep("HK-AML", times = length(hkaml_sum)), rep("Beat-AML", times = length(beataml_sum)), rep("TCGA-LAML", times = length(tcga_sum))), count = c(hkaml_sum, beataml_sum, tcga_sum), stringsAsFactors = FALSE)
df_tmm$dataset <- as.factor(df_tmm$dataset)
df_tmm$dataset <- factor(df_tmm$dataset, levels = c("HK-AML", "TCGA-LAML", "Beat-AML"))

# Figure 1a
ggplot(df_tmm, aes(x=dataset, y=count, fill = dataset)) +
   geom_boxplot() +
   theme_minimal() +
   xlab("AML dataset") +
   ylab("TMM at SE loci") +
   theme01

# Figure 1b
hkaml_nonZero <- colSums(hkaml_count>0)
beataml_nonZero <- colSums(beataml_count>0)
tcga_nonZero <- colSums(tcga_count>0)

number_sEnh_df <- data.frame(sEnh_number = c(hkaml_nonZero, tcga_nonZero, beataml_nonZero),
                             dataset = c(rep("HK-AML", times=190), rep("TCGA-LAML", times=144), rep("Beat-AML", times=510)))
number_sEnh_df$dataset <- factor(number_sEnh_df$dataset, levels = c("HK-AML", "TCGA-LAML", "Beat-AML"))

ggplot(number_sEnh_df, aes(x=dataset, y=sEnh_number, fill = dataset)) +
   geom_boxplot() +
   theme_minimal() +
   xlab("AML dataset") +
   ylab("Number of transcribed SE loci") +
   ylim(0,300000) +
   theme01

rm(hkaml_count, beataml_count, tcga_count)
rm(hkaml_lib.size, beataml_lib.size, tcga_lib.size)
rm(hkaml_sum, beataml_sum, tcga_sum)
rm(hmrf_nonZero, tcga_nonZero, beataml_nonZero)
rm(dge.hkaml, dge.beataml, dge.tcga)

# Figure 1c
# Read batch-corrected transcribed SE matrix for HK-AML cohort. This matrix includes only the 77563 significantly expressed SE loci.

df_globalEnhancerActivation <- data.frame(superEnh = colSums(combat_superEnh))
df_globalEnhancerActivation$Type <- as.factor(c(rep("AML", times=185), rep("Control", times=5)))

ggplot(data = df_globalEnhancerActivation, aes(x=Type, y=superEnh, fill=Type)) +
   geom_boxplot() +
   xlab("Specimen Type") +
   ylab("Expression of SE Loci") +
   ylim(0,100000) +
   theme_minimal() +
   geom_signif(annotations = paste0("p=", formatC(wilcox.test(colSums(combat_superEnh[,1:185]), colSums(combat_superEnh[,186:190]))$p.value, digits = 2)),
               y_position = 95000, xmin=1, xmax=2, tip_length = c(500000, 500000), textsize = 8) +
   theme01

rm(df_globalEnhancerActivation)

# Figure 1d
# Group SE loci into SE regions (Hnisz et al. Cell 2013) for hierarchical clustering
download.file("https://figshare.com/ndownloader/files/43118452?private_link=081f4f46d84c084b0098", "./tmp_files/gr.loci_sEnh.rds")
gr.loci_sEnh <- readRDS("./tmp_files/gr.loci_sEnh.rds")

hnisz_rownames <- vector()
hnisz_sum_logTMM <- list()
hnisz_sEnh_length <- vector()
j <- 1
for (hniszLoci in unique(gr.loci_sEnh$sEnh_region_name)){
   sEnh_temp <- gr.loci_sEnh$sEnh[which(gr.loci_sEnh$sEnh_region_name == hniszLoci)]
   sEnh_temp <- sEnh_temp[sEnh_temp %in% rownames(combat_superEnh)]
   if (length(sEnh_temp) > 1){
      df_temp <- combat_superEnh[rownames(combat_superEnh) %in% sEnh_temp,]
      hnisz_sum_logTMM[[j]] <- apply(df_temp, 2, sum)
      hnisz_rownames[j] <- hniszLoci
      hnisz_sEnh_length[j] <- length(sEnh_temp)
      j <- j + 1
   }
   else if (length(sEnh_temp) == 1) {
      df_temp <- combat_superEnh[rownames(combat_superEnh) %in% sEnh_temp,]
      hnisz_sum_logTMM[[j]] <- df_temp
      hnisz_rownames[j] <- hniszLoci
      hnisz_sEnh_length[j] <- length(sEnh_temp)
      j <- j + 1
   }
   else{
      # Do nothing  
   }
}
rm(df_temp)
df_hnisz_sum_logTMM <- as.data.frame(matrix(unlist(hnisz_sum_logTMM), nrow = length(hnisz_sum_logTMM), byrow = TRUE))
rownames(df_hnisz_sum_logTMM) <- hnisz_rownames
colnames(df_hnisz_sum_logTMM) <- colnames(combat_superEnh)
pairwise_cor_hnisz_bySample_allHnisz <- rcorr(as.matrix(df_hnisz_sum_logTMM), type = "spearman")

# Download clinical classification and mutational data for heatmap annotation

idx_rga <- group_final$Group == "t(8;21)" | group_final$Group == "inv(16)/t(16;16)" | group_final$Group == "KMT2A" | group_final$Group == "bZIP CEBPA" | group_final$Group == "NPM1"
idx_rga <- c(idx_rga, rep(FALSE, times = 5))

group_rga <- group_final[idx_rga, c(15,16,18)]
group_rga$Group <- as.factor(as.character(group_rga$Group))

annotation_colour01 <- brewer.pal(n = 9, "Paired")
names(annotation_colour01) <- levels(group_final$Group)
annotation_colour02 <- c("#D3D3D3", "#4682B4")
names(annotation_colour02) <- levels(group_final$DNMT3A)

annotation_colour_all <- list(
   Group = annotation_colour01,
   DNMT3A = annotation_colour02,
   `FLT3 ITD` = annotation_colour02,
   ASXL1 = annotation_colour02,
   EZH2 = annotation_colour02,
   BCOR = annotation_colour02,
   STAG2 = annotation_colour02,
   SRSF2 = annotation_colour02,
   SF3B1 = annotation_colour02,
   U2AF1 = annotation_colour02,
   ZRSR2 = annotation_colour02,
   RUNX1 = annotation_colour02,
   `Aneuploidy/TP53` = annotation_colour02,
   TET2 = annotation_colour02,
   IDH1 = annotation_colour02,
   `IDH2 p.R140` = annotation_colour02,
   `IDH2 p.R172` = annotation_colour02
)
annotation_colour_reorder <- annotation_colour_all
annotation_colour_reorder$WHO <- annotation_colour_reorder$Group[-c(3)]
names(annotation_colour_reorder$WHO)[1] <- "AML-MR"

annotation_colour_rga <- list(
   Group = annotation_colour01,
   DNMT3A = annotation_colour02,
   `FLT3 ITD` = annotation_colour02
)
annotation_colour_rga$Group <- annotation_colour_rga$Group[c(2,4,5,6,9)]

pheatmap(pairwise_cor_hnisz_bySample_allHnisz$r[idx_rga,idx_rga],
         annotation_col = group_rga[,1:3],
         clustering_method = "average",
         annotation_colors = annotation_colour_rga,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = F,
         show_colnames = F,
         color = rev(colorRampPalette(brewer.pal(n=7, name = "RdBu"))(100)),
         breaks = seq(0.7, 1.0, by=0.003))

# Figure S1a
pheatmap(pairwise_cor_hnisz_bySample_allHnisz$r[1:185,1:185],
         annotation_col = group_final,
         clustering_method = "average",
         annotation_colors = annotation_colour_reorder,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = F,
         show_colnames = F,
         color = rev(colorRampPalette(brewer.pal(n=7, name = "RdBu"))(100)),
         breaks = seq(0.7, 1.0, by=0.003))

rm(annotation_colour_all, annotation_colour_reorder, annotation_colour_rga, annotation_colour01, annotation_colour02)
rm(group_rga, pairwise_cor_hnisz_bySample_allHnisz, hnisz_rownames, hnisz_sEnh_length, hniszLoci, j)

# Figure 1e
# Plot the relationship between Spearman rho and distance between TSS and SE loci (stratified into within-loop and cross-boundary in TAD)


ggplot(cor_sEnh_ENSG_2Mb, aes(x = dist_TSS_sEnh, y = rho, colour = ENSG_in_TAD)) +
   geom_smooth() +
   xlim(0,1000) +
   theme_minimal() +
   xlab("Distance (kb)") +
   ylab("Correlation (Spearman rho)") +
   scale_color_discrete(name = "SE Locus-Gene Pairs", labels = c("Cross TAD boundary", "Within TAD loop")) +
   theme(axis.text = element_text(size = 22), axis.title = element_text(size = 24), legend.title = element_text(size = 22), legend.text = element_text(size = 20))

# Figure S1b

sEnh_intergenic <- geneRegion_sEnh$sEnh[geneRegion_sEnh$variable == "intergenic"]
cor_sEnh_ENSG_intergenic <- cor_sEnh_ENSG_2Mb[cor_sEnh_ENSG_2Mb$sEnh %in% sEnh_intergenic,]

ggplot(cor_sEnh_ENSG_intergenic, aes(x = dist_TSS_sEnh, y = rho, colour = ENSG_in_TAD)) +
   geom_smooth() +
   xlim(0,1000) +
   theme_minimal() +
   xlab("Distance (kb)") +
   ylab("Correlation (Spearman rho)") +
   scale_color_discrete(name = "SE Locus-Gene Pairs", labels = c("Cross TAD boundary", "Within TAD loop")) +
   theme(axis.text = element_text(size = 22), axis.title = element_text(size = 24), legend.title = element_text(size = 22), legend.text = element_text(size = 20))

rm(cor_sEnh_ENSG_final, cor_sEnh_ENSG_2Mb, cor_sEnh_ENSG_intergenic)

# Figure 1f

enrich_AML_RGA <- enrich_AML_RGA[grepl("(AML|ACUTE_MYELOID_LEUKEMIA|MULLIGHAN)", rownames(enrich_AML_RGA)),]
enrich_AML_RGA[is.na(enrich_AML_RGA)] <- 1   # Replace NA with 1
enrich_AML_RGA <- -log10(enrich_AML_RGA)

# Rearrange the rows according to AML subtypes
idx <- c(44,52,74,42,45,62,73,28,30,27,29,46,48,54,63,72,75,2,1,32,34,36,31,33,35,68,67,57,64)
enrich_AML_RGA <- enrich_AML_RGA[idx,]
colnames(enrich_AML_RGA) <- c("t(8;21)", "inv(16)", "KMT2A", "NPM1", "CEBPA")
rownames(enrich_AML_RGA) <- gsub("_", " ", rownames(enrich_AML_RGA))

pheatmap(enrich_AML_RGA[c(18,20:22,26,5:7,4,1:3,8:9,12,13:16,28,29),],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         angle_col = 45,
         color = rev(colorRampPalette(brewer.pal(n=7, name = "RdBu"))(100)),
         fontsize = 14
)