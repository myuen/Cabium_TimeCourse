library(Mfuzz)

library(dplyr)
library(stringr)
library(tibble)
library(tidyselect)

# Function to calculate row variance
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

cpm <- read.table("results/tmm_normalized.cpm.txt")

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
    
    # Calculate mean for treated and untreated libraries
    "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
    "T.mean" = rowMeans(select(cpm, starts_with("T"))),
    
    # Calculate variance for treated and untreated libraries
    "MJ.var" = rowVar(select(cpm, starts_with("M"))),
    "T.var" = rowVar(select(cpm, starts_with("T")))
)


# Take only the most highly express / variable
# This is just for the sake of this tutorial
# test_data <-
#   cpm.stats[which(cpm.stats$MJ.var > 50 &
#                     cpm.stats$MJ.mean > 50), c(1:5)]
# str(test_data)
# 'data.frame':	3997 obs. of  7 variables:

mj.top <- cpm.stats %>%
  filter(MJ.mean > 100 & MJ.var > 100) %>%
  select(1:5)

mj.top <- cpm.stats %>%
  select(1:5)

mj.top <- mj.top[!rowSums(mj.top) == 0,]

str(mj.top)

# all.equal(test_data, mj.top)
# [1] TRUE


time <- c(2, 24, 48, 96, 192)

mj.top <- rbind(time, mj.top)

rownames(mj.top)[1] <- "time"

tmp <- tempfile()

write.table(mj.top, file = tmp, sep = '\t', quote = F, col.names = NA)

data <- table2eset(tmp)

data.s <- standardise(data)

m1 <- mestimate(data.s)
# [1] 2.037463

Dmin(data.s, m = m1, crange = seq(2, 22, 1), repeats = 5, visu = TRUE)
# Crashed with full set of data

clust = 8

c <- mfuzz(data.s, c = clust, m = m1)

# mfuzz.plot(data.s, cl = c, mfrow = c(4, 3),
#            time.labels = c(2, 24, 48, 96, 192),
#            new.window = FALSE)


mfuzz.plot2(data.s, cl = c, mfrow = c(4, 3), colo = "fancy",
           time.labels = c(2, 24, 48, 96, 192),
           new.window = FALSE)

acore <- acore(data.s, c, min.acore = 0)

acore_list <- 
  do.call(rbind, lapply(seq_along(acore), function(i){
    data.frame(CLUSTER=i, acore[[i]])}))


### Clustering base on logFC instead of CPM (?)

