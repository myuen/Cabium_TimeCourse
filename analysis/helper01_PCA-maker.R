# argument 1, expression matrix
# argument 2, groupings, as factor


PCA_maker <- function(cpm, groupings) {
  require(ggplot2)
  
  pca <- prcomp(t(cpm), scale. = TRUE)
  pca_summary <- summary(pca)
  pc <- as.data.frame(pca$x)
  
  g <- guide_legend("Genotype")
  
  p <-
    ggplot(
      data = pc,
      aes_string(x = "PC1", y = "PC2", shape = factor(groupings))) +
    geom_point(size = 4) +
    
    # Naming title and axis
    ggtitle("Principal Component Analysis") + 
    xlab(paste("PC1 (", pca_summary$importance["Proportion of Variance", "PC1"]*100, "%)")) +
    ylab(paste("PC2 (", pca_summary$importance["Proportion of Variance", "PC2"]*100, "%)")) +
    
    # Text on labels
    geom_text(data = pc, aes(x = PC1, y = PC2, label = rownames(pc)),
              size = 3, vjust = 2.25, hjust = 0.5) +
    
    theme_bw() + 
    
    # Setting legend colour and shape
    guides(colour = guide_legend("Genotype"), shape = guide_legend("Genotype")) +
    
    scale_colour_manual(values = c("firebrick1", "dodgerblue3"))
  
  return(p)
}
