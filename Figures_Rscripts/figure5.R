# R codes for generating Figure 5 in the transcribed SE paper
library(ggplot2)
library(Seurat)

options(timeout = 60000)

# Figure 5a

ggplot(geneRegion_sEnh_2data, aes(x = variable, fill=dataNature)) +
   geom_bar(position = "dodge") +
   coord_flip() +
   ylab("Count") +
   xlab("Gene Region") +
   scale_x_discrete(limits = rev(levels(geneRegion_sEnh_2data$variable))) +
   theme_minimal() +
   scale_fill_discrete(name = "Data Type") +
   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 20), legend.text = element_text(size = 20))

# Figure 5b
# Download the result files from Seurat. Each around 1 Gb.

# Plot the annotation of cell types in NPM1
DimPlot(AML556all, label = TRUE, repel = TRUE, label.size = 6) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle("NPM1-mutated") +
        theme(plot.title = element_text(hjust = 0.5))

# Plot the pooled SE loci expression by cells
FeaturePlot(AML556all, features = "nCount_SCT_sEnh", slot = 'scale.data', max.cutoff = 200, cols = c("steelblue", "red")) +
   ggtitle("NPM1-mutated") +
   xlab("UMAP 1") +
   ylab("UMAP 2")

# Figure 5c
# Plot the annotation of cell types in t(8;21)
DimPlot(AML707Ball, label = TRUE, repel = TRUE, label.size = 6) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        ggtitle("t(8;21)") +
        theme(plot.title = element_text(hjust = 0.5))

# Plot the pooled SE loci expression by cells
FeaturePlot(AML707Ball, features = "nCount_SCT_sEnh", slot = 'scale.data', max.cutoff = 50, cols = c("steelblue", "red")) +
   ggtitle("t(8;21)") +
   xlab("UMAP 1") +
   ylab("UMAP 2")
