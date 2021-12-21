library(edgeR)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tximport)


source("analysis/99_helper01_PCA-maker.R")


### Differential Expression Analysis on tween (control) and methyl-jasmonate
### (treated) laser-micro-dissected cambium celltype over 5 time points
### in PG-653 white spruce

# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


# Import file with tximport
quant.files <-
  dir("data",
      pattern = "quant.sf",
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE)

samples <-
  str_extract(quant.files, "[MT]\\S+rep\\d")

names(quant.files) <- samples


# Read Salmon quantification with tximport
txi <-
  tximport(quant.files,
           type = "salmon",
           txOut = TRUE,
           varReduce = TRUE,
           countsFromAbundance = "lengthScaledTPM")

treatment <- str_split(samples, "-")


# Create experimental design
expDes <- 
  data.frame(
    sample = samples,
    treatment = factor(samples %>% str_extract("[M|T]\\w"),
                       levels = c("TW", "MJ")),
    timepoint = as.character(samples %>% 
                         str_extract("T\\d") %>% 
                         str_replace("T", "")),
    biorep = as.character(samples %>% str_extract("\\d$")))

expDes$group <- with(expDes, interaction(treatment, timepoint))

#        sample treatment timepoint biorep group
# 1  MJ-T1_rep1        MJ         1      1  MJ.1
# 2  MJ-T1_rep2        MJ         1      2  MJ.1
# 3  MJ-T1_rep3        MJ         1      3  MJ.1
# 4  MJ-T2_rep1        MJ         2      1  MJ.2
# 5  MJ-T2_rep2        MJ         2      2  MJ.2
# 6  MJ-T2_rep3        MJ         2      3  MJ.2
# 7  MJ-T3_rep1        MJ         3      1  MJ.3
# 8  MJ-T3_rep2        MJ         3      2  MJ.3
# 9  MJ-T3_rep3        MJ         3      3  MJ.3
# 10 MJ-T4_rep1        MJ         4      1  MJ.4
# 11 MJ-T4_rep2        MJ         4      2  MJ.4
# 12 MJ-T4_rep3        MJ         4      3  MJ.4
# 13 MJ-T5_rep1        MJ         5      1  MJ.5
# 14 MJ-T5_rep2        MJ         5      2  MJ.5
# 15 MJ-T5_rep3        MJ         5      3  MJ.5
# 16 TW-T1_rep1        TW         1      1  TW.1
# 17 TW-T1_rep2        TW         1      2  TW.1
# 18 TW-T1_rep3        TW         1      3  TW.1
# 19 TW-T2_rep1        TW         2      1  TW.2
# 20 TW-T2_rep2        TW         2      2  TW.2
# 21 TW-T2_rep3        TW         2      3  TW.2
# 22 TW-T3_rep1        TW         3      1  TW.3
# 23 TW-T3_rep2        TW         3      2  TW.3
# 24 TW-T3_rep3        TW         3      3  TW.3
# 25 TW-T4_rep1        TW         4      1  TW.4
# 26 TW-T4_rep2        TW         4      2  TW.4
# 27 TW-T4_rep3        TW         4      3  TW.4
# 28 TW-T5_rep1        TW         5      1  TW.5
# 29 TW-T5_rep2        TW         5      2  TW.5
# 30 TW-T5_rep3        TW         5      3  TW.5


# Load counts into DGEList object from edgeR package.
x <- DGEList(counts = txi$counts, 
             group = expDes$group)

dim(x)
# [1] 47338    30


# Filter low-expression contigs.  Keep only genes with at 
# least 1 count-per-million reads (cpm) in at least 3 samples
# (i.e. at least expressed in 1 timepoint in all biological replicates)
y <- x[(rowSums(cpm(x) > 1) >= 3), ]

dim(y)
# [1] 30129    30


# Reset depth
y$samples$lib.size <- colSums(y$counts)


# TMM Normalization by Depth
y <- calcNormFactors(y)


# make model matrix
# Nested design model matrix
modMat <- model.matrix(~treatment + treatment:timepoint, expDes)

colnames(modMat) <- colnames(modMat) %>% 
  str_replace(":", "_") %>% 
  str_replace("treatment", "") %>% 
  str_replace("timepoint", "T")

#    (Intercept) MJ TW_T2 MJ_T2 TW_T3 MJ_T3 TW_T4 MJ_T4 TW_T5 MJ_T5
# 1            1  1     0     0     0     0     0     0     0     0
# 2            1  1     0     0     0     0     0     0     0     0
# 3            1  1     0     0     0     0     0     0     0     0
# 4            1  1     0     1     0     0     0     0     0     0
# 5            1  1     0     1     0     0     0     0     0     0
# 6            1  1     0     1     0     0     0     0     0     0
# 7            1  1     0     0     0     1     0     0     0     0
# 8            1  1     0     0     0     1     0     0     0     0
# 9            1  1     0     0     0     1     0     0     0     0
# 10           1  1     0     0     0     0     0     1     0     0
# 11           1  1     0     0     0     0     0     1     0     0
# 12           1  1     0     0     0     0     0     1     0     0
# 13           1  1     0     0     0     0     0     0     0     1
# 14           1  1     0     0     0     0     0     0     0     1
# 15           1  1     0     0     0     0     0     0     0     1
# 16           1  0     0     0     0     0     0     0     0     0
# 17           1  0     0     0     0     0     0     0     0     0
# 18           1  0     0     0     0     0     0     0     0     0
# 19           1  0     1     0     0     0     0     0     0     0
# 20           1  0     1     0     0     0     0     0     0     0
# 21           1  0     1     0     0     0     0     0     0     0
# 22           1  0     0     0     1     0     0     0     0     0
# 23           1  0     0     0     1     0     0     0     0     0
# 24           1  0     0     0     1     0     0     0     0     0
# 25           1  0     0     0     0     0     1     0     0     0
# 26           1  0     0     0     0     0     1     0     0     0
# 27           1  0     0     0     0     0     1     0     0     0
# 28           1  0     0     0     0     0     0     0     1     0
# 29           1  0     0     0     0     0     0     0     1     0
# 30           1  0     0     0     0     0     0     0     1     0


# voom transformation
v <- voom(y, modMat, plot = TRUE)

# Transformed expression value
v.cpm <- v$E

# Write voom transformed expression matrix
write.table(
  v.cpm,
  "results/cpm.txt",
  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Create PCA plot
# p <- PCA_maker(v.cpm, 
#                colnames(v.cpm) %>% str_replace("_rep\\d", ""))

# ggsave("results/figures/PCA.jpg", plot = p, dpi = "retina")


# Linear modeling
fit <- lmFit(v, modMat)

cont_matrix <-
  makeContrasts(
    MJ.T1_TW.T1 = MJ,
    MJ.T2_TW.T2 = MJ_T2 - TW_T2,
    MJ.T3_TW.T3 = MJ_T3 - TW_T3,
    MJ.T4_TW.T4 = MJ_T4 - TW_T4,
    MJ.T5_TW.T5 = MJ_T5 - TW_T5,
    levels = modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

# Up-regulation = higher expression in gland in reference to glandes
summary(
  decideTests(fit3, method = "separate", 
              adjust.method = "fdr", p.value = pCutoff, 
              lfc = lfcCutoff))

#        MJ.T1_TW.T1 MJ.T2_TW.T2 MJ.T3_TW.T3 MJ.T4_TW.T4 MJ.T5_TW.T5
# Down           278        1039        2008        1584        1657
# NotSig       29346       28151       26805       27349       27020
# Up             505         939        1316        1196        1452

# results <-
#   topTable(fit2, coef = "gland", 
#            sort.by = "logFC", number = Inf)

# str(results)


# Reorganize columns
# results <- results %>% 
#   select(logFC, adj.P.Val) %>%
#   rownames_to_column("cds")


### Write out stats of all contigs
# write.table(
#   results, quote = FALSE, sep = "\t", 
#   row.names = FALSE, col.names = TRUE,
#   "results/dea_all_stats.txt"
# )

### Write out all sig DE stats
# sigDE <- results %>% 
#   filter(adj.P.Val <= pCutoff & abs(logFC) >= lfcCutoff)
# 
# write.table(
#   sigDE, quote = FALSE, sep = "\t", 
#   row.names = FALSE, col.names = TRUE,
#   "results/dea_sigDE_stats.txt"
# )
# 
