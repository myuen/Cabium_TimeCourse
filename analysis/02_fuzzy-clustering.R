library(Mfuzz)

library(dplyr)
library(tibble)


# Read normalized CPM
cpm <- read.table("results/tmm_normalized.cpm.txt")


# Function to calculate row variance
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


cpm.stats <-
  cpm %>% transmute(
    "MJ.T1" = rowMeans(select(cpm, starts_with("MJ.T1"))),
    "MJ.T2" = rowMeans(select(cpm, starts_with("MJ.T2"))),
    "MJ.T3" = rowMeans(select(cpm, starts_with("MJ.T3"))),
    "MJ.T4" = rowMeans(select(cpm, starts_with("MJ.T4"))),
    "MJ.T5" = rowMeans(select(cpm, starts_with("MJ.T5"))),
    
    "TW.T1" = rowMeans(select(cpm, starts_with("TW.T1"))),
    "TW.T2" = rowMeans(select(cpm, starts_with("TW.T2"))),
    "TW.T3" = rowMeans(select(cpm, starts_with("TW.T3"))),
    "TW.T4" = rowMeans(select(cpm, starts_with("TW.T4"))),
    "TW.T5" = rowMeans(select(cpm, starts_with("TW.T5"))),
    
    "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
    "TW.mean" = rowMeans(select(cpm, starts_with("T"))),
    
    "MJ.var" = rowVar(select(cpm, starts_with("M"))),
    "TW.var" = rowVar(select(cpm, starts_with("T")))
  )


# Read statistics for differential expression sequences
sigDE <- read.table("results/dea_sigDE_stats.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE)

str(sigDE)
# 'data.frame':	12915 obs. of  4 variables:


# Filter CPM if it is sig. differentially expressed
sigDE.stats <- 
  cpm.stats[rownames(cpm.stats) %in% sigDE$cds,] %>% 
  filter(MJ.var != 0) %>% 
  filter(TW.var != 0) %>% 
  select(1:5)

str(sigDE.stats)
# 'data.frame':	6762 obs. of  5 variables:


# Add time to table
time <- c(2, 24, 48, 96, 192)

sigDE.stats <- rbind(time, sigDE.stats)

rownames(sigDE.stats)[1] <- "time"


# table2eset function must read from a file
write.table(sigDE.stats, file = "results/sigDE.stats.tbl.txt",
            sep = '\t', quote = F, col.names = NA)


# Convert an expression table to eset object
data <- table2eset("results/sigDE.stats.tbl.feb17.txt")


# Standardise data
data.s <- standardise(data)


(m1 <- mestimate(data.s))
# [1] 2.005024


pdf("results/figures/dmin.pdf", width = 20, height = 10)

Dmin(data.s, m = m1, crange = seq(2, 10, 1), repeats = 5, visu = TRUE)
# [1] 2.2974304 2.3323899 2.0215294 1.1607365 1.1254128 1.0166034 0.7747944
# [8] 0.8162585 0.7447335

dev.off()


# Number of soft clusters
clust <- 5

c <- mfuzz(data.s, c = clust, m = m1)

###

clusters <- c[[1]]
clusterindex <- c[[3]]
memship <- c[[4]]
memship[memship < 0.7] <- -1
colorindex <- integer(dim(exprs(data.s))[[1]])

x <- matrix(rnorm(5*5),nrow=5) # new expression matrix with two genes
mem.tmp <- membership(x, clusters=clusters,m=ml) #membership values

###


# Write out standardized expression data
s.data <- exprs(data.s) %>% 
  as.data.frame() %>% 
  rownames_to_column('cds')

write.table(s.data, "results/standardized_expression.txt", 
            quote = FALSE, row.names = FALSE)
            
# Write out cluster members
clstr_members <- as.data.frame(c[3])

clstr_members <- clstr_members %>% 
  rownames_to_column("cds")

write.table(clstr_members, "results/cluster_members.feb17.txt", 
            quote = FALSE, row.names = FALSE)
