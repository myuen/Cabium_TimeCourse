library(Mfuzz)

library(dplyr)
library(stringr)
library(tibble)
library(tidyselect)

<<<<<<< HEAD
cpm <- read.table("results/tmm_normalized.cpm.txt")

# cpm <- read.table("~/TPM_Pg653_Cambium_TRD_4april2017_nr.txt")

rowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

=======
# Function to calculate row variance
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

cpm <- read.table("results/tmm_normalized.cpm.txt")

>>>>>>> 59ea49d104c283cfeb52ad68283c58f7732c97f3
cpm.stats <-
  cpm %>% transmute(
    "MJ.T1" = rowMeans(select(cpm, starts_with("MJ.T1"))),
    "MJ.T2" = rowMeans(select(cpm, starts_with("MJ.T2"))),
    "MJ.T3" = rowMeans(select(cpm, starts_with("MJ.T3"))),
    "MJ.T4" = rowMeans(select(cpm, starts_with("MJ.T4"))),
    "MJ.T5" = rowMeans(select(cpm, starts_with("MJ.T5"))),
<<<<<<< HEAD

=======
>>>>>>> 59ea49d104c283cfeb52ad68283c58f7732c97f3
    "TW.T1" = rowMeans(select(cpm, starts_with("TW.T1"))),
    "TW.T2" = rowMeans(select(cpm, starts_with("TW.T2"))),
    "TW.T3" = rowMeans(select(cpm, starts_with("TW.T3"))),
    "TW.T4" = rowMeans(select(cpm, starts_with("TW.T4"))),
    "TW.T5" = rowMeans(select(cpm, starts_with("TW.T5"))),
<<<<<<< HEAD

    "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
    "TW.mean" = rowMeans(select(cpm, starts_with("T"))),

    "MJ.var" = rowVar(select(cpm, starts_with("M"))),
    "T.var" = rowVar(select(cpm, starts_with("T")))
  )


# cpm.stats <-
#   cpm %>% transmute(
#     "MJ.T1" = rowMeans(select(cpm, starts_with("M") & ends_with("2H"))),
#     "MJ.T2" = rowMeans(select(cpm, starts_with("M") & ends_with("1D"))),
#     "MJ.T3" = rowMeans(select(cpm, starts_with("M") & ends_with("2D"))),
#     "MJ.T4" = rowMeans(select(cpm, starts_with("M") & ends_with("4D"))),
#     "MJ.T5" = rowMeans(select(cpm, starts_with("M") & ends_with("8D"))),
# 
#     "TW.T1" = rowMeans(select(cpm, starts_with("T") & ends_with("2H"))),
#     "TW.T2" = rowMeans(select(cpm, starts_with("T") & ends_with("1D"))),
#     "TW.T3" = rowMeans(select(cpm, starts_with("T") & ends_with("2D"))),
#     "TW.T4" = rowMeans(select(cpm, starts_with("T") & ends_with("4D"))),
#     "TW.T5" = rowMeans(select(cpm, starts_with("T") & ends_with("8D"))),
# 
#     "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
#     "TW.mean" = rowMeans(select(cpm, starts_with("T"))),
# 
#     "MJ.var" = rowVar(select(cpm, starts_with("M"))),
#     "T.var" = rowVar(select(cpm, starts_with("T")))
#     )

sigDE <- read.table("results/dea_sigDE_stats.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	12915 obs. of  4 variables:


#take only the most highly express/ variable
#this is just for the sake of this tutorial
# mj.top50 <- cpm.stats %>% 
#   filter(MJ.mean > 50 & MJ.var > 50 ) %>%
#   select(1:5)
# 
# str(mj.top50)
# 'data.frame':	2705 obs. of  5 variables:

# T2 <- sigDE %>% filter(focus == "MJT2vsTWT2")
# 
# T2.stats <- cpm.stats[rownames(cpm.stats) %in% T2$cds,] %>% 
#   filter(MJ.var != 0) %>% 
#   select(1:5)

# str(T2.stats)
# 'data.frame':	2390 obs. of  5 variables:

sigDE.stats <- 
  cpm.stats[rownames(cpm.stats) %in% sigDE$cds,] %>% 
  filter(MJ.var != 0) %>% 
  select(1:5)

str(sigDE.stats)
# 'data.frame':	6960 obs. of  5 variables:

time <- c(2, 24, 48, 96, 192)

sigDE <- rbind(time, sigDE)

rownames(sigDE)[1] <- "time"

tmp <- tempfile()

write.table(sigDE, file = tmp, sep = '\t', quote = F, col.names = NA)
=======
    
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
>>>>>>> 59ea49d104c283cfeb52ad68283c58f7732c97f3

data <- table2eset(tmp)

data.s <- standardise(data)

<<<<<<< HEAD
(m1 <- mestimate(data.s))
# [1] 2.004545

Dmin(data.s, m = m1, crange = seq(2, 10, 1), repeats = 5, visu = TRUE)

clust <- 5

c <- mfuzz(data.s, c = clust, m = m1)

mfuzz.plot(
  data.s,
  cl = c,
  mfrow = c(2, 4),
  time.labels = c(2, 24, 48, 96, 192),
  new.window = FALSE
)


svg("results/5-clusters-for-dea.svg", width = 20, height = 10)

mfuzz.plot2(data.s, cl = c, mfrow = c(2, 4), colo = "fancy",
            time.labels = c(2, 24, 48, 96, 192),
            xlab = "Time Point (hrs)",
            centre = TRUE, 
            x11 = FALSE)
dev.off()

clusters <- as.data.frame(c[3])
clusters <- clusters %>% rownames_to_column("cds")
=======
m1 <- mestimate(data.s)
# [1] 2.037463

Dmin(data.s, m = m1, crange = seq(2, 22, 1), repeats = 5, visu = TRUE)
# Crashed with full set of data

clust = 8

c <- mfuzz(data.s, c = clust, m = m1)

# mfuzz.plot(data.s, cl = c, mfrow = c(4, 3),
#            time.labels = c(2, 24, 48, 96, 192),
#            new.window = FALSE)

>>>>>>> 59ea49d104c283cfeb52ad68283c58f7732c97f3

write.table(clusters, "results/cluster_members.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

acore <- acore(data.s, c, min.acore = 0)

acore_list <- 
  do.call(rbind, lapply(seq_along(acore), function(i){
    data.frame(CLUSTER=i, acore[[i]])}))


### Clustering base on logFC instead of CPM (?)

