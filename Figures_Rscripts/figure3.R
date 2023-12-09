# R codes for generating Figures 3d to 3f in the transcribed SE paper
# Figures 3a to 3b are tracks, while Figure 3c is an output diagram from dmrseq, and are not included here
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
options(timeout = 6000)
options(stringsAsFactors = FALSE)

# Define general theme for plotting
theme03a <- theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 18),
                 legend.title = element_text(size = 20), legend.text = element_text(size = 18), legend.position = "none")
theme03b <- theme(axis.text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size = 20), plot.title = element_text(size = 18), legend.position = "none", plot.margin = unit(c(1,1,1,3), "cm"))

# Figure 3d

# Plot scatter plot for the 4 SEs in the vicinity of BCL6
p <- list()
i <- 1
for (k in c(21168, 21179, 21186, 21198)){
   temp_hnisz <- paste0("hnisz", k)
   df_temp <- data.frame(hnisz = df_hnisz_median_logTMM[temp_hnisz, 1:185], BCL6 = exp_main[1:185,"ENSG00000113916"], DNMT3A = as.factor(group_final$DNMT3A[1:185]))
   p[[i]] <- ggplot(df_temp, aes(y = BCL6, x = hnisz)) +
      geom_point(aes(colour = DNMT3A), cex = 2) +
      ylab("BCL6 Expression") +
      xlab(paste0("Expression of ", temp_hnisz)) +
      stat_cor(method = "pearson", label.x = 0.2, label.y = 9.5, size = 6) +
      theme_minimal() +
      ylim(0, 10) +
      theme03a
   i <- i+1
}
do.call(gridExtra::grid.arrange, arrangeGrob(grobs = p, ncol = 4))
rm(df_temp, p, i, k)

# Figure 3e
# Work on an dataframe for data visualization
temp_df_BCL6 <- data.frame(Exp = exp_main[,"ENSG00000113916"], hnisz21198 = df_hnisz_median_logTMM["hnisz21198",], group = group_final$Group, DNMT3A = group_final$DNMT3A)
temp_df_BCL6 <- temp_df_BCL6[temp_df_BCL6$group != "Control",]
temp_df_BCL6$final_group <- as.character(temp_df_BCL6$group)
temp_df_BCL6$final_group[temp_df_BCL6$group == "NPM1" & temp_df_BCL6$DNMT3A == "Positive"] <- "NPM1+ DNMT3A+"
temp_df_BCL6$final_group[temp_df_BCL6$group == "NPM1" & temp_df_BCL6$DNMT3A == "Negative"] <- "NPM1+ DNMT3A-"
temp_df_BCL6$final_group[temp_df_BCL6$group == "Biallelic_CEBPA"] <- "CEBPA"
temp_df_BCL6$final_group <- as.factor(temp_df_BCL6$final_group)
temp_df_BCL6$final_group <- factor(temp_df_BCL6$final_group, levels = levels(temp_df_BCL6$final_group)[c(7,6,4,10,1,8,9,3,5,2)])

# Plot barplot on SE expression (hnisz21198) in AML subtypes
ggplot(temp_df_BCL6, aes(x=final_group, y=hnisz21198, fill=final_group)) +
   geom_boxplot(show.legend = FALSE, colour = "blue", fill="royalblue", alpha=0.2) +
   theme_minimal() +
   xlab("AML Subtype") +
   ylab("SE (hnisz21198) Expression") +
   theme03b

# Figure 3f
# Plot barplot on BCL6 expression in AML subtypes
ggplot(temp_df_BCL6, aes(x=final_group, y=Exp)) +
   geom_boxplot(show.legend = FALSE, colour = "blue", fill="royalblue", alpha=0.2) +
   theme_minimal() +
   xlab("AML Subtype") +
   ylab("BCL6 Expression") +
   ylim(0,8) +
   theme03b

rm(temp_df_BCL6)