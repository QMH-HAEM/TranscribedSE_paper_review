# R codes for generating Figure 6 in the transcribed SE paper
library(pheatmap)
library(RColorBrewer)

# Download the short listed drug response data

# Plot heatmap
pheatmap(drugCor_hnisz_shortlist,
         labels_col = temp_labels,
         fontsize_col = 10,
         fontsize_row = 12,
         angle_col = 45,
         color = rev(colorRampPalette(brewer.pal(n=7, name = "RdBu"))(100))
         )
