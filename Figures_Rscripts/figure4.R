# R codes for generating Figures 4a to 4b in the transcribed SE paper
# Figure 4c is a diagram on tracks, and is not included here
library(dplyr)
library(grid)
library(pheatmap)
library(RColorBrewer)

options(timeout = 6000)
options(stringsAsFactors = FALSE)

# Figure 4a

# Group genes with same interaction(s) with SEs
temp_label_hnisz_DE_heatmap <- as.data.frame(df_hnisz_pcg_final2[,c(6,3)] %>% distinct())
label_hnisz_DE_heatmapFinal2 <- list()
j <- 1
for (i in 1:nrow(temp_label_hnisz_DE_heatmap)){
   if (i == 1){
      label_hnisz_DE_heatmapFinal2[[j]] <- temp_label_hnisz_DE_heatmap[i,]
   }
   else{
      if (temp_label_hnisz_DE_heatmap[i,1] == label_hnisz_DE_heatmapFinal2[[j-1]][1]){
         j <- j - 1
         label_hnisz_DE_heatmapFinal2[[j]][2] <- paste0(label_hnisz_DE_heatmapFinal2[[j]][2], ",", temp_label_hnisz_DE_heatmap[i,2])
      }
      else{
         label_hnisz_DE_heatmapFinal2[[j]] <- temp_label_hnisz_DE_heatmap[i,]
      }
   }
   j <- j + 1
}
label_hnisz_DE_heatmapFinal2 <- data.frame(hnisz = sapply(label_hnisz_DE_heatmapFinal2, '[[', 1), label = sapply(label_hnisz_DE_heatmapFinal2, '[[', 2))
temp_label_PCG_heatmap <- label_hnisz_DE_heatmapFinal2

# Shorten some labels of genes for easier visualization
label_hnisz_DE_heatmapFinal2[255,2] <- "HOXA genes,SNX10"
label_hnisz_DE_heatmapFinal2[256,2] <- "HOXA genes,HOTAIRM1"
label_hnisz_DE_heatmapFinal2[257,2] <- "HOXA genes,HOTTIP"
label_hnisz_DE_heatmapFinal2[819,2] <- "HBG1,HBG2,HBD,HBB"
label_hnisz_DE_heatmapFinal2[179,2] <- "NFIA"
label_hnisz_DE_heatmapFinal2[1128,2] <- "HMGA2"
label_hnisz_DE_heatmapFinal2[3011,2] <- "PRDM16"
label_hnisz_DE_heatmapFinal2[3056,2] <- "PRDM16"
label_hnisz_DE_heatmapFinal2[706,2] <- "MPP7,MPP7-DT"
label_hnisz_DE_heatmapFinal2[538,2] <- "MAMDC2"
label_hnisz_DE_heatmapFinal2[738,2] <- "LRMDA,ZNF503"
label_hnisz_DE_heatmapFinal2[2394,2] <- "LTBP1"
label_hnisz_DE_heatmapFinal2[2898,2] <- "HHIP"
label_hnisz_DE_heatmapFinal2[2513,2] <- "SLC38A11,GRB14,COBLL1,SCN3A"
label_hnisz_DE_heatmapFinal2[1148,2] <- "SOCS2,CRADD"
label_hnisz_DE_heatmapFinal2[2707,2] <- "HTR1F"
label_hnisz_DE_heatmapFinal2[2851,2] <- "MMRN1,SNCA"
label_hnisz_DE_heatmapFinal2[545,2] <- "SLC25A24,FAM102B"
label_hnisz_DE_heatmapFinal2[2917,2] <- "TMA16,TKTL2"
label_hnisz_DE_heatmapFinal2[166,2] <- "TACSTD2"
label_hnisz_DE_heatmapFinal2[2101,2] <- "SAMD11"
label_hnisz_DE_heatmapFinal2[2567,2] <- "IKZF2"
label_hnisz_DE_heatmapFinal2[2568,2] <- "IKZF2"
label_hnisz_DE_heatmapFinal2[2569,2] <- "IKZF2"
label_hnisz_DE_heatmapFinal2[2570,2] <- "IKZF2"
label_hnisz_DE_heatmapFinal2[671,2] <- "SFMBT2"
label_hnisz_DE_heatmapFinal2[1250,2] <- "SMIM2,TUSC8"
label_hnisz_DE_heatmapFinal2[2202,2] <- "APP"
label_hnisz_DE_heatmapFinal2[2203,2] <- "APP"
label_hnisz_DE_heatmapFinal2[1808,2] <- "CD300C/E/LD/LF"
label_hnisz_DE_heatmapFinal2[1836,2] <- "BAHCC1"
label_hnisz_DE_heatmapFinal2[2436,2] <- "MEIS1"
label_hnisz_DE_heatmapFinal2[2437,2] <- "MEIS1"

temp_label <- label_hnisz_DE_heatmapFinal2$label[match(rownames(df_hnisz_median_logTMM), label_hnisz_DE_heatmapFinal2$hnisz)]
names(temp_label) <- label_hnisz_DE_heatmapFinal2$hnisz[match(rownames(df_hnisz_median_logTMM), label_hnisz_DE_heatmapFinal2$hnisz)]
temp_label[is.na(temp_label)] <- ""

# Give annotation bar as in Figure S1a
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

# Plot heatmap
pheatmap(df_hnisz_median_logTMM[rownames(df_hnisz_median_logTMM) %in% subtypeSpec_hnisz, cluster_order],
         annotation_col = group_final[,17, drop = FALSE],
         annotation_colors = annotation_colour_reorder,
         cluster_cols = FALSE,
         clustering_method = "complete",
         labels_row = temp_label[rownames(df_hnisz_median_logTMM) %in% subtypeSpec_hnisz],
         show_colnames = F,
         scale = "row",
         color = rev(colorRampPalette(brewer.pal(n=7, name = "RdBu"))(100)),
         breaks = seq(-6, 6, by=0.12)
)

# Figure 4b
# Load the result tables merged from matchIt returns

pval2_adj.peaks_haemCells <- matrix(p.adjust(pval2.peaks_haemCells, method = "fdr"), nrow = 8, byrow = F)
colnames(pval2_adj.peaks_haemCells) <- colnames(pval2.peaks_haemCells)
rownames(pval2_adj.peaks_haemCells) <- rownames(pval2.peaks_haemCells)
colnames(odds2.peaks_haemCells) <- gsub("_", "-", colnames(odds2.peaks_haemCells))
rownames(odds2.peaks_haemCells) <- c("AML", "t(8;21)", "inv(16)", "KMT2A", "NPM1", "CEBPA", "AML-MR", "NPM1cluster1vs23")

# Adjust labeling of X-axis
draw_colnames_45 <- function (coln, gaps, ...) {
   coord <- pheatmap:::find_coordinates(length(coln), gaps)
   x     <- coord$coord - 0.5 * coord$size
   res   <- grid::textGrob(
      coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
      vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
   )
   return(res)
}
assignInNamespace(
   x = "draw_colnames",
   value = "draw_colnames_45",
   ns = asNamespace("pheatmap")
)

# Plot heatmap for the heptad factors
pheatmap(t(odds2.peaks_haemCells[2:7,c(5:24,33:36,41:44)]),
         cluster_rows = F,
         cluster_cols = T
)
