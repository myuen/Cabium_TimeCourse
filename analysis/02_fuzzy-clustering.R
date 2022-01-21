library(Mfuzz)

library(dplyr)
library(stringr)
library(tibble)
library(tidyselect)


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


sigDE.stats <- 
  cpm.stats[rownames(cpm.stats) %in% sigDE$cds,] %>% 
  filter(MJ.var != 0) %>% 
  select(1:5)

str(sigDE.stats)
# 'data.frame':	6960 obs. of  5 variables:

time <- c(2, 24, 48, 96, 192)

sigDE.stats <- rbind(time, sigDE.stats)

rownames(sigDE.stats)[1] <- "time"


# table2eset function must read from a file
write.table(sigDE.stats, file = "results/sigDE.stats.tbl.txt", 
            sep = '\t', quote = F, col.names = NA)


data <- table2eset("results/sigDE.stats.tbl.txt")

data.s <- standardise(data)

(m1 <- mestimate(data.s))
# [1] 2.004545

pdf("results/dmin.pdf", width = 20, height = 10)

tmp <- Dmin(data.s, m = m1, crange = seq(2, 10, 1), repeats = 5, visu = TRUE)
# [1] 2.2943019 2.3378306 2.0339792 1.0750610 1.1141968 1.0213772 0.9464997 0.8278212 0.7318490

dev.off()


clust <- 5

c <- mfuzz(data.s, c = clust, m = m1)

pdf("results/five-clusters-for-dea.pdf", width = 20, height = 10)

mfuzz.plot2(data.s, cl = c, mfrow = c(2, 3), colo = "fancy",
            time.labels = c(2, 24, 48, 96, 192),
            xlab = "Time (hrs)",
            ylab = expression("Log"[2]*'FC'),
            centre = TRUE, 
            # Increase font size of title
            cex.main = 1.2,
            # Increase font size of tick point
            cex.axis = 1.1,
            # Increase font size of axis lable
            cex.lab = 1.2,
            x11 = FALSE)

dev.off()


# Write out cluster members
clusters <- as.data.frame(c[3])

clusters <- clusters %>% rownames_to_column("cds")

write.table(clusters, "results/cluster_members.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
