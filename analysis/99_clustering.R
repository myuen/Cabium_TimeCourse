library(Mfuzz)
library(Biobase)

library(dplyr)
library(stringr)
library(tibble)
library(tidyselect)


# cpm <- read.table("results/tmm_normalized.cpm.txt")

# mean.cpm <- 
#   cpm %>% transmute(
#     "MJ.T1" = rowMeans(select(cpm, starts_with("MJ.T1"))),
#     "MJ.T2" = rowMeans(select(cpm, starts_with("MJ.T2"))),
#     "MJ.T3" = rowMeans(select(cpm, starts_with("MJ.T3"))),
#     "MJ.T4" = rowMeans(select(cpm, starts_with("MJ.T4"))),
#     "MJ.T5" = rowMeans(select(cpm, starts_with("MJ.T5"))),
#     "TW.T1" = rowMeans(select(cpm, starts_with("TW.T1"))),
#     "TW.T2" = rowMeans(select(cpm, starts_with("TW.T2"))),
#     "TW.T3" = rowMeans(select(cpm, starts_with("TW.T3"))),
#     "TW.T4" = rowMeans(select(cpm, starts_with("TW.T4"))),
#     "TW.T5" = rowMeans(select(cpm, starts_with("TW.T5"))))

rowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

cpm <- read.table("~/TPM_Pg653_Cambium_TRD_4april2017_nr.txt")

cpm.stats <- 
  cpm %>% transmute(
    "MJ.T1" = rowMeans(select(cpm, starts_with("M") & ends_with("2H"))),
    "MJ.T2" = rowMeans(select(cpm, starts_with("M") & ends_with("1D"))),
    "MJ.T3" = rowMeans(select(cpm, starts_with("M") & ends_with("2D"))),
    "MJ.T4" = rowMeans(select(cpm, starts_with("M") & ends_with("4D"))),
    "MJ.T5" = rowMeans(select(cpm, starts_with("M") & ends_with("8D"))),
    
    "T.T1" = rowMeans(select(cpm, starts_with("T") & ends_with("2H"))),
    "T.T2" = rowMeans(select(cpm, starts_with("T") & ends_with("1D"))),
    "T.T3" = rowMeans(select(cpm, starts_with("T") & ends_with("2D"))),
    "T.T4" = rowMeans(select(cpm, starts_with("T") & ends_with("4D"))),
    "T.T5" = rowMeans(select(cpm, starts_with("T") & ends_with("8D"))),
    
    "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
    "T.mean" = rowMeans(select(cpm, starts_with("T"))),
    
    "MJ.var" = rowVar(select(cpm.stats, starts_with("M"))),
    "T.var" = rowVar(select(cpm.stats, starts_with("T")))
    )


#take only the most highly express/ variable
#this is just for the sake of this tutorial
mj.top50 <- cpm.stats %>% 
  filter(MJ.mean > 50 & MJ.var > 50) %>% 
  select(1:5)
str(mj.top50)
# 'data.frame':	2707 obs. of  5 variables:

time <- c(2, 24, 48, 96, 192)

mj.top50 <- rbind(time, mj.top50)

rownames(mj.top50)[1] <- "time"

tmp <- tempfile()

write.table(mj.top50, file = tmp, sep = '\t', quote = F, col.names = NA)

data <- table2eset(tmp)

data.s <- standardise(data)

m1 <- mestimate(data.s)
# [1] 2.026014

Dmin(data.s, m=m1, crange = seq(2, 22, 1), repeats = 3, visu = TRUE)

#  [1] 2.1524100 1.6868709 1.7339129 1.8031554 1.3578340 0.9946767 1.0649987 0.9930474 0.9366941 0.7816127
# [11] 0.6810765 0.6781544 0.6765621 0.6571941 0.6441000 0.6193729 0.5877256 0.5739728 0.5480030 0.5602557
# [21] 0.5320178

clust = 12

c <- mfuzz(data.s, c = clust, m = m1)

mfuzz.plot(data.s, cl = c, mfrow = c(1, 1),
           time.labels = c(2, 24, 48, 96, 192),
           new.window = FALSE)

mfuzz.plot2(data.s, cl = c, mfrow = c(4, 3), colo = "fancy",
           time.labels = c(2, 24, 48, 96, 192),
           new.window = FALSE)

