# R codes for generating Figure 2 in the transcribed SE paper
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(ggVennDiagram)
library(org.Hs.eg.db)
options(timeout = 6000)
options(stringsAsFactors = FALSE)
options(scipen = 3)

# Define general theme for plotting
theme02 <- theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 18), legend.position = "none")

# Figure 2a

# Group sources of interaction data
df_hnisz_pcg_final2$Source <- "SE database"
df_hnisz_pcg_final2$Source[df_hnisz_pcg_final2$dataOrigin=="assi_NPM1" |
                              df_hnisz_pcg_final2$dataOrigin=="assi_t821" |
                              df_hnisz_pcg_final2$dataOrigin=="ptasinska" |
                              df_hnisz_pcg_final2$dataOrigin=="mifsud"|
                              df_hnisz_pcg_final2$dataOrigin=="hic_CEBPA" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_FLT3" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_HSPC" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_inv16" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_KMT2B" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_NPM1" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_RUNX1" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_STAG2" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_t69" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_t821" |
                              df_hnisz_pcg_final2$dataOrigin=="hic_TET2" |
                              df_hnisz_pcg_final2$dataOrigin=="hichip_CMP" |
                              df_hnisz_pcg_final2$dataOrigin=="hichip_GMP" |
                              df_hnisz_pcg_final2$dataOrigin=="hichip_HSC" |
                              df_hnisz_pcg_final2$dataOrigin=="hichip_MEP"] <- "Hi-C based"
df_hnisz_pcg_final2$Source[df_hnisz_pcg_final2$dataOrigin=="EBI" | df_hnisz_pcg_final2$dataOrigin=="GTEx" | df_hnisz_pcg_final2$dataOrigin=="pancanqtl"] <- "eQTL"

df_hnisz_pcg_final2$idx <- paste0(df_hnisz_pcg_final2$sEnh, "_", df_hnisz_pcg_final2$ENSG)

# Plot Venn diagram
venn_input2 <- list(
   `SE database` = df_hnisz_pcg_final2$idx[df_hnisz_pcg_final2$Source=="SE database"],
   `Hi-C` = df_hnisz_pcg_final2$idx[df_hnisz_pcg_final2$Source=="Hi-C based"],
   `eQTL` = df_hnisz_pcg_final2$idx[df_hnisz_pcg_final2$Source=="eQTL"]
)
ggVennDiagram(venn_input2) +
   scale_color_manual(values = c("black", "black", "black"))
#   scale_fill_distiller(palette = palette(c("white", "white", "white")))

venn <- Venn(venn_input2)
data <- process_data(venn)
ggplot() +
   geom_sf(aes(color = id, lwd = 6), data = venn_setedge(data), show.legend = FALSE) +
   geom_sf_text(aes(label = name, size = 36), data = venn_setlabel(data)) +
   geom_sf_label(aes(label = count, size = 36), data = venn_region(data)) +
   theme_void() +
   theme(legend.position = "none")
rm(venn, venn_input2, data)

# Figure 2b
# Load SE loci annotation and protein annotation from Gencode

gr.anno_ENSG <- sortSeqlevels(gr.anno_ENSG)
gr.anno_ENSG <- sort(gr.anno_ENSG)

gr.loci_sEnh$nearest_ENSG_row <- nearest(gr.loci_sEnh, gr.anno_ENSG)
gr.loci_sEnh$nearest_gene_name <- gr.anno_ENSG$external_gene_name[nearest(gr.loci_sEnh, gr.anno_ENSG)]
gr.loci_sEnh$nearest_ENSG <- gr.anno_ENSG$ensembl_gene_id[nearest(gr.loci_sEnh, gr.anno_ENSG)]

anno_sEnh_forClosestNeighbour <- anno_sEnh_cor[,1:2] %>% unique()   # 57566

# If ENSG exists in gr.anno_ENSG, return row number in gr.anno_ENSG
# Compare returned row number with gr.loci_sEnh$nearest_ENSG_row
# Return difference in row number in a list for each sEnh
# Sort the list and obtain the smallest difference
# If zero, closest neighbouring protein-coding gene. Otherwise not.
gr.loci_sEnh$returnFromRowComparison <- list(c(99999))
for (i in 1:nrow(anno_sEnh_forClosestNeighbour)){
   if (anno_sEnh_forClosestNeighbour$ENSG[i] %in% gr.anno_ENSG$ensembl_gene_id){
      temp_row_num <- which(gr.anno_ENSG$ensembl_gene_id == anno_sEnh_forClosestNeighbour$ENSG[i])
      temp_diff <- temp_row_num - gr.loci_sEnh$nearest_ENSG_row[gr.loci_sEnh$sEnh == anno_sEnh_forClosestNeighbour$sEnh[i]]
      gr.loci_sEnh$returnFromRowComparison[gr.loci_sEnh$sEnh == anno_sEnh_forClosestNeighbour$sEnh[i]] <- list(c(unlist(gr.loci_sEnh$returnFromRowComparison[gr.loci_sEnh$sEnh == anno_sEnh_forClosestNeighbour$sEnh[i]]), temp_diff))
   }
   else {
      next
   }
}

hist_data <- table(sapply(sapply(gr.loci_sEnh$returnFromRowComparison, abs), min))
hist_data <- hist_data[-41]   # Exclude the entries with 99999, i.e. no match
hist_data[4] <- sum(hist_data[4:40])
hist_data <- hist_data[1:4]
names(hist_data)[1] <- "Closest"
names(hist_data)[4] <- ">= 3"

df_temp <- data.frame(`Position` = hist_data)
colnames(df_temp) <- c("Gene", "Position")
#df_temp 

# Plot distribution
ggplot(data=df_temp, aes(x=Gene, y=Position)) +
   geom_bar(stat="identity", fill = "steelblue") +
   theme_minimal() +
   theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), plot.title = element_text(size = 20), legend.position = "none") +
   xlab("Position of Interacting Gene") +
   ylab("Frequency")

# Figure 2c: Volcano plot for differentially expressed SEs in NPM1-mutated AML

download.file("https://figshare.com/ndownloader/files/43148872?private_link=7502326da13d9f870a77", "./tmp_files/gencode30_proteinCoding_list.rds")
gencode30_proteinCoding_list <- readRDS("./tmp_files/gencode30_proteinCoding_list.rds")

# Work on labeling of SEs using their linked protein-coding genes
df_hnisz_pcg_plot$genelabels_NPM1 <- FALSE
df_hnisz_pcg_plot$genelabels_NPM1[grepl("^HOXA\\d+$", df_hnisz_pcg_plot$gene)] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "MEIS1"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "BAHCC1"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "PRDM16"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "IKZF2"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "CD300E"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "LTBP1"] <- TRUE
df_hnisz_pcg_plot$genelabels_NPM1[df_hnisz_pcg_plot$gene == "CD34" & df_hnisz_pcg_plot$DE_hnisz_NPM1_fdr < 0.01 & (df_hnisz_pcg_plot$DE_hnisz_NPM1_logFC > 0.6 | df_hnisz_pcg_plot$DE_hnisz_NPM1_logFC < -0.6)] <- TRUE

temp <- df_hnisz_pcg_plot %>%
   as_tibble() %>%
   filter(ENSG %in% gencode30_proteinCoding_list) %>%
   mutate(DE = ((DE_hnisz_NPM1_fdr < 0.01)&
                   DE_sEnh_NPM1_fdr < 0.05 &
                   DE_ENSG_NPM1_fdr < 0.05 &
                   (DE_hnisz_NPM1_logFC > 2 | DE_hnisz_NPM1_logFC < -1.6) &
                   (sign(DE_ENSG_NPM1_logFC) == sign(DE_hnisz_NPM1_logFC)))
   ) %>%
   mutate(DE2 = (DE_hnisz_NPM1_fdr < 0.01 & (DE_hnisz_NPM1_logFC > 0.6 | DE_hnisz_NPM1_logFC < -0.6))) %>%
   dplyr::select(ENSG, gene, hnisz, DE_hnisz_NPM1_fdr, DE_hnisz_NPM1_logFC, DE_ENSG_NPM1_fdr, DE_ENSG_NPM1_logFC, genelabels_NPM1, DE, DE2) %>%
   unique()

temp$label <- NA
temp$label[temp$genelabels_NPM1 & temp$DE] <- temp$gene[temp$genelabels_NPM1 & temp$DE]

# Plot volcano plot
ggplot(data = temp, aes(x=DE_hnisz_NPM1_logFC, y=-log10(DE_hnisz_NPM1_fdr), col=DE2, label=label)) +
   geom_point(size = 3) +
   xlab("log2 fold change") +
   ylab("-log10 FDR") +
   xlim(-12, 12) +
   theme_minimal() +
   geom_text_repel(colour = "black", box.padding = 0.5, max.overlaps = Inf, size = 6) +
   theme02
rm(temp)

# Figure 2d: Volcano plot for differentially expressed SEs in KMT2A-rearranged AML
# Work on labeling of SEs using their linked protein-coding genes
df_hnisz_pcg_plot$genelabels_KMT2A <- FALSE
df_hnisz_pcg_plot$genelabels_KMT2A[grepl("^HOXA\\d+$", df_hnisz_pcg_plot$gene)] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "MEIS1"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "MEF2C"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "MLLT10"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "MECOM"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "ZNF521"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "COL5A1"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "TRPS1"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "SNCA"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "EREG"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "MMRN1"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "PLIN2"] <- TRUE
df_hnisz_pcg_plot$genelabels_KMT2A[df_hnisz_pcg_plot$gene == "HTR1F"] <- TRUE

temp <- df_hnisz_pcg_plot %>%
   as_tibble() %>%
   filter(ENSG %in% gencode30_proteinCoding_list) %>%
   mutate(DE = ((DE_hnisz_KMT2A_fdr < 0.001 | (((DE_hnisz_KMT2A_logFC > 3 | DE_hnisz_KMT2A_logFC < -3) & DE_hnisz_KMT2A_fdr < 0.01))) &
                   DE_sEnh_KMT2A_fdr < 0.05 &
                   DE_ENSG_KMT2A_fdr < 0.05 &
                   (DE_hnisz_KMT2A_logFC > 0.6 | DE_hnisz_KMT2A_logFC < -0.6) &
                   (sign(DE_ENSG_KMT2A_logFC) == sign(DE_hnisz_KMT2A_logFC)))
   ) %>%
   mutate(DE2 = (DE_hnisz_KMT2A_fdr < 0.01 & (DE_hnisz_KMT2A_logFC > 0.6 | DE_hnisz_KMT2A_logFC < -0.6))) %>%
   dplyr::select(ENSG, gene, hnisz, DE_hnisz_KMT2A_fdr, DE_hnisz_KMT2A_logFC, DE_ENSG_KMT2A_fdr, DE_ENSG_KMT2A_logFC, genelabels_KMT2A, DE, DE2) %>%
   unique()

temp$label <- NA
temp$label[temp$genelabels_KMT2A & temp$DE] <- temp$gene[temp$genelabels_KMT2A & temp$DE]

# Plot volcano plot
ggplot(data = temp, aes(x=DE_hnisz_KMT2A_logFC, y=-log10(DE_hnisz_KMT2A_fdr), col=DE2, label=label)) +
   geom_point(size = 3) +
   xlab("log2 fold change") +
   ylab("-log10 FDR") +
   xlim(-12, 12) +
   ylim(0,6.5) +
   theme_minimal() +
   geom_text_repel(colour = "black", box.padding = 0.5, max.overlaps = Inf, size = 6) +
   theme02
rm(temp)

# Figure 2e
# Load the lists of protein-coding genes associated with differentially expressed SE loci in AML subtypes

# Define function to extract Entrez ID from Ensembl gene ID for downstream over-representation analysis
ENSGtoEntrez2 <- function(ENSG){
   entrezList <- mapIds(org.Hs.eg.db,
                        keys = ENSG, 
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")
   return(entrezList)
}

# Define an empty master Entrez ID list to fill in
df_entrez_AML <- data.frame(Entrez=character(0), Expression=character(0), Subtype=character(0))

# t(8;21)
temp_entrez_hyper <- ENSGtoEntrez2(pcg_t821_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "t(8;21)")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# inv(16)
temp_entrez_hyper <- ENSGtoEntrez2(pcg_inv16_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "inv(16)")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# KMT2A
temp_entrez_hyper <- ENSGtoEntrez2(pcg_KMT2A_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "KMT2A")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# NPM1
temp_entrez_hyper <- ENSGtoEntrez2(pcg_NPM1_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "NPM1")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# CEBPA
temp_entrez_hyper <- ENSGtoEntrez2(pcg_CEBPA_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "CEBPA")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# TP53
temp_entrez_hyper <- ENSGtoEntrez2(pcg_TP53_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "TP53")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# Chromatin-spliceosome
temp_entrez_hyper <- ENSGtoEntrez2(pcg_ChromSplice_sEnhHyper$ENSG) 
df_entrez_temp <- data.frame(Entrez=c(temp_entrez_hyper),
                             Expression=c(rep("Overexpressed", times = length(temp_entrez_hyper))),
                             Subtype = "ChrSpli")
df_entrez_AML <- rbind(df_entrez_AML, df_entrez_temp)
rm(temp_entrez_hyper, df_entrez_temp)

# Perform over-representation analysis and visualisation
ora_compare_AML_GO <- compareCluster(Entrez~Expression+Subtype, data = df_entrez_AML, OrgDb = org.Hs.eg.db, fun = "enrichGO", ont = "BP")
clusterProfiler::dotplot(ora_compare_AML_GO, showCategory=5, x=~Expression) + facet_grid(~Subtype) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Figure 2f
# Load the matched data

# Over-expressed SE loci & chromatin pattern in DNMT3A-mutated AML
atac_table <- table(m.data.randomMatchedPairs$atac_DNMT3A_specific, m.data.randomMatchedPairs$hnisz_hyper)
atac_fisher <- fisher.test(atac_table)

temp_df <- data.frame(
   SE = factor(c("Over-expressed", "Not over-expressed"), levels = c("Over-expressed", "Not over-expressed")),
   Percentage = c(atac_table[2,2] / sum(atac_table[,2]) * 100,
                  atac_table[2,1] / sum(atac_table[,1]) * 100)
)
ggplot(temp_df, aes(x=SE, y=Percentage)) +
   geom_bar(stat = "identity", fill = "steelblue") +
   xlab("Super-enhancer") +
   ylab("Overlap open chromatin (%)") +
   ylim(0,110) +
   theme_minimal() +
   theme(axis.text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size = 18), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
   geom_signif(annotations = paste0("p=", atac_fisher$p.value), y_position = 100, xmin=1, xmax=2, tip_length = c(5, 5), textsize = 8)
rm(temp_df, atac_table)

# Figure 2g
# Over-expressed SE loci & DNA methylation in DNMT3A-mutated AML
dmr_table <- table(m.data.randomMatchedPairs$dmrseq_hypo, m.data.randomMatchedPairs$hnisz_hyper)
dmr_fisher <- fisher.test(dmr_table)

temp_df <- data.frame(
   SE = factor(c("Over-expressed", "Not over-expressed"), levels = c("Over-expressed", "Not over-expressed")),
   Percentage = c(dmr_table[2,2] / sum(dmr_table[,2]) * 100,
                  dmr_table[2,1] / sum(dmr_table[,1]) * 100)
)
ggplot(temp_df, aes(x=SE, y=Percentage)) +
   geom_bar(stat = "identity", fill = "steelblue") +
   xlab("Super-enhancer") +
   ylab("Overlap hypo-DMR (%)") +
   ylim(0,45) +
   theme_minimal() +
   theme(axis.text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size = 18), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
   geom_signif(annotations = paste0("p=", dmr_fisher$p.value), y_position = 41, xmin=1, xmax=2, tip_length = c(0.2, 0.2), textsize = 8)
rm(temp_df, dmr_table)

# Figure 2h
df_hnisz_pcg_plot$genelabels_NPM1cluster1vs23[df_hnisz_pcg_plot$DE_hnisz_NPM1cluster1vs23_fdr < 0.005] <- TRUE

temp <- df_hnisz_pcg_plot %>%
   as_tibble() %>%
   filter(ENSG %in% gencode30_proteinCoding_list) %>%
   mutate(DE = (DE_hnisz_NPM1cluster1vs23_fdr < 0.01 &
                   DE_sEnh_NPM1cluster1vs23_fdr < 0.05 &
                   DE_ENSG_NPM1cluster1vs23_fdr < 0.05 &
                   DE_hnisz_NPM1cluster1vs23_logFC > 0.6 &
                   (sign(DE_ENSG_NPM1cluster1vs23_logFC) == sign(DE_hnisz_NPM1cluster1vs23_logFC)))
   ) %>%
   mutate(DE2 = (DE_hnisz_NPM1cluster1vs23_fdr < 0.01 & (DE_hnisz_NPM1cluster1vs23_logFC > 0.6 | DE_hnisz_NPM1cluster1vs23_logFC < -0.6))) %>%
   dplyr::select(ENSG, gene, hnisz, DE_hnisz_NPM1cluster1vs23_fdr, DE_hnisz_NPM1cluster1vs23_logFC, DE_ENSG_NPM1cluster1vs23_fdr, DE_ENSG_NPM1cluster1vs23_logFC, genelabels_NPM1cluster1vs23, DE, DE2) %>%
   unique()

temp$label <- NA
temp$label[temp$genelabels_NPM1cluster1vs23 & temp$DE] <- temp$gene[temp$genelabels_NPM1cluster1vs23 & temp$DE]
temp$label[!(temp$label %in% gobp_adaptive_immune_response)] <- NA
temp$label[temp$DE_hnisz_NPM1cluster1vs23_logFC > 2 & temp$DE_hnisz_NPM1cluster1vs23_fdr < 0.000001 & temp$DE] <- temp$gene[temp$DE_hnisz_NPM1cluster1vs23_logFC > 2 & temp$DE_hnisz_NPM1cluster1vs23_fdr < 0.000001 & temp$DE]
temp_CD300 <- which(grepl("CD300", temp$label))
temp <- temp[-temp_CD300[-1],]
temp$label[temp_CD300[1]] <- "CD300A/C/E/H/LB/LD/LF"

ggplot(data = temp, aes(x=DE_hnisz_NPM1cluster1vs23_logFC, y=-log10(DE_hnisz_NPM1cluster1vs23_fdr), col=DE2, label=label)) +
   geom_point(size = 3) +
   xlab("log2 fold change") +
   ylab("-log10 FDR") +
   xlim(-12, 12) +
   theme_minimal() +
   geom_text_repel(colour = "black", box.padding = 0.5, max.overlaps = Inf, size = 6) +
   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 18), legend.position = "none")

# Figure 2i
pcg_NPM1cluster1vs23_hniszHyper <- df_hnisz_pcg_final2 %>%
   as_tibble() %>%
   filter(ENSG %in% gencode30_proteinCoding_list) %>%
   filter(DE_sEnh_NPM1cluster1vs23_fdr < 0.01) %>%
   filter(DE_sEnh_NPM1cluster1vs23_logFC > 0) %>%
   filter(DE_ENSG_NPM1cluster1vs23_fdr < 0.05) %>%
   filter(DE_ENSG_NPM1cluster1vs23_logFC > 0) %>%
   dplyr::select(ENSG, gene, hnisz, DE_hnisz_NPM1cluster1vs23_fdr, DE_hnisz_NPM1cluster1vs23_logFC, DE_ENSG_NPM1cluster1vs23_fdr, DE_ENSG_NPM1cluster1vs23_logFC) %>%
   unique()

# Download GO gene sets (obtained from MSigDB)

   # Convert ENSG number to Entrez ID
   anno_genelist <- mapIds(org.Hs.eg.db,
                           keys = genelist, 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
   universe_genelist <- mapIds(org.Hs.eg.db,
                               keys = universe, 
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")
   
   # Perform enrichment study using msigdb
   msigdb <- read.gmt(gmt)
   enrich <- enricher(as.character(anno_genelist), universe = as.character(universe_genelist), TERM2GENE = msigdb, qvalueCutoff = 0.05)
   return(enrich)
}
enrichGO_NPM1cluster1vs23_hniszHyper <- msigDB_enrich_2(genelist = unique(pcg_NPM1cluster1vs23_hniszHyper$ENSG), universe = unique(gencode30_proteinCoding_list))
dotplot(enrichGO_NPM1cluster1vs23_hniszHyper, showCategory = 9, font.size = 14)
