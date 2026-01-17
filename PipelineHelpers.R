nicePCA <- function(pcaOut, groupLogic) {
  # PCA data
  pca_df <- as.data.frame(pcaOut$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  
  # logical group to labeled factor
  pca_df$Group <- factor(ifelse(groupLogic, "2", "1"))
  
  # variance explained %
  pve <- round(pcaOut$pve[1:2] * 100, 2)
  
  # point plot
  plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw(base_size = 14) +
    labs(
      x = paste0("PCA 1 (", pve[1], "% variance explained)"),
      y = paste0("PCA 2 (", pve[2], "% variance explained)")
    ) +
    coord_equal() 
  
  return(plot)
}


niceVolcano <- function(resOutput, pAdjCutoff, log2FCcutoff = 1.5) {
  # count for subtitle in plot later
  lociCount <- nrow(resOutput)
  
  # remove rows w/o padj
  resOutput <- resOutput[!is.na(resOutput$padj), ]
  
  # points that exceed both thresholds
  resOutput$ExceedsBoth <- factor((resOutput$padj < pAdjCutoff) & (abs(resOutput$log2FoldChange) > log2FCcutoff))
  
  # significance col
  resOutput$SigDE <- factor(resOutput$padj < pAdjCutoff)
  
  # top 5 up and down genes exceeding both cutoffs
  topUp <- head(resOutput[resOutput$ExceedsBoth == TRUE & resOutput$log2FoldChange > 0, ][
    order(resOutput$padj[resOutput$ExceedsBoth == TRUE & resOutput$log2FoldChange > 0]), ], 5)
  topDown <- head(resOutput[resOutput$ExceedsBoth == TRUE & resOutput$log2FoldChange < 0, ][
    order(resOutput$padj[resOutput$ExceedsBoth == TRUE & resOutput$log2FoldChange < 0]), ], 5)
  labelGenes <- rbind(topUp, topDown)
  
  # volcano plot
  p <- ggplot(resOutput, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = ExceedsBoth, shape = SigDE), alpha = 1, size = 2) +
    
    # threshold lines
    geom_vline(xintercept = -log2FCcutoff, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = log2FCcutoff, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(pAdjCutoff), linetype = "dashed", color = "black") +
    
    # labs for top genes w/ lines
    geom_text_repel(
      data = labelGenes,
      aes(label = gene_name),
      max.overlaps = 20,
      segment.color = "black",
      segment.size = 0.5,
      min.segment.length = 0.1,
      nudge_y = 0.2
    ) +
    
    # labs
    labs(
      x = "log2(Fold Change)",
      y = "−log10(FDR−Adjusted p−values)",
      subtitle = paste(lociCount, "loci analyzed")
    ) +
    
    # Reorder legends: SigDE (shape) first, ExceedsBoth (color) second
    guides(
      shape = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      plot.subtitle = element_text(hjust = 0, vjust = 1)
    )
  
  return(p)
}


resSummary <- function(resOutput) {

  totalGenes <- nrow(resOutput)
  print("How many genes were included in the DESeq2 analysis after filtering?")
  print(totalGenes)
  
  nonMissingPadj <- sum(!is.na(resOutput$padj))
  print("How many genes included in the DESeq2 analysis had non-missing p-adjusted values?")
  print(nonMissingPadj)
  
  print("Note: using 0.05 as the p-adjusted cutoff for significance for the following calculations.")
  sigCutoff <- 0.05
  
  downregulated <- sum(!is.na(resOutput$padj) & resOutput$padj < sigCutoff & resOutput$log2FoldChange < 0)
  print("How many genes were significantly downregulated in group2 relative to group1?")
  print(downregulated)
  
  upregulated <- sum(!is.na(resOutput$padj) & resOutput$padj < sigCutoff & resOutput$log2FoldChange > 0)
  print("How many genes were significantly upregulated in group2 relative to group1?")
  print(upregulated)
  
  topGenes <- head(resOutput[order(resOutput$padj), ], 10)
  topGeneNames <- topGenes$gene_name
  print("Names of the top 10 most significantly differentially expressed genes from smallest to largest p-adjusted:")
  print(topGeneNames)
}
