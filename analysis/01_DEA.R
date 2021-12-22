library(edgeR)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tximport)


source("analysis/helper01_PCA-maker.R")


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
  dir("data/Salmon",
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
x <- DGEList(counts = txi$counts, group = expDes$group)

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


# With nested design model matrix, we can test the effect of 
# MJ treatment on specific timepoint
modMat <- model.matrix(~treatment * timepoint, expDes)

colnames(modMat)[1] <- "Intercept"

colnames(modMat) <- colnames(modMat) %>% str_replace(":", "_")


# Inspecting the differential expression from treatment at 
# different time-point without the compounding effect
contrasts <-
  makeContrasts(
    MJT1vsTWT1 = treatmentMJ,
    # The effect of MJ treatment at T2 without the time effect
    MJT2vsTWT2 = treatmentMJ + treatmentMJ_timepoint2,
    # The effect of MJ treatment at T3 without the time effect
    MJT3vsTWT3 = treatmentMJ + treatmentMJ_timepoint3,
    # The effect of MJ treatment at T4 without the time effect
    MJT4vsTWT4 = treatmentMJ + treatmentMJ_timepoint4,
    # The effect of MJ treatment at T5 without the time effect
    MJT5vsTWT5 = treatmentMJ + treatmentMJ_timepoint5,
    levels = modMat)

# voom transformation
v <- voom(y, modMat, plot = TRUE)

# Transformed expression value
v.cpm <- v$E

# Write voom transformed expression matrix
# write.table(
#   v.cpm, "results/cpm.txt",
#   row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Create PCA plot
# p <- PCA_maker(v.cpm,
#                colnames(v.cpm) %>% str_replace("_rep\\d", ""))

# ggsave("results/figures/PCA.jpg", plot = p, dpi = "retina")


# Linear modeling
fit <- lmFit(v, modMat)

fit2 <- contrasts.fit(fit, contrasts)

fit3 <- eBayes(fit2)

# Up-regulation = higher expression in gland in reference to glandes
summary(
  decideTests(fit3, method = "separate",
              adjust.method = "fdr", p.value = pCutoff, 
              lfc = lfcCutoff))
#        MJT1vsTWT1 MJT2vsTWT2 MJT3vsTWT3 MJT4vsTWT4 MJT5vsTWT5
# Down          278       1018       1749       1530       1309
# NotSig      29346      27738      26568      26869      27209
# Up            505       1373       1812       1730       1611


# Get the results with topTable for each coef
rslt <- map_dfr(colnames(contrasts), function(x) {
  df <- topTable(fit3, coef = x, number = Inf, sort.by = "none")
  df <- df %>% rownames_to_column("cds")
  df <- df %>% select(cds, logFC, adj.P.Val)
  df <- df %>% add_column("focus" = x)
})

table(rslt$focus)
# MJT1vsTWT1 MJT2vsTWT2 MJT3vsTWT3 MJT4vsTWT4 MJT5vsTWT5 
#      30129      30129      30129      30129      30129 

# Write out stats of all contigs
write.table(
  rslt, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/dea_all_stats.txt"
)

# Write out all sig DE stats
sigDE <- rslt %>%
  filter(adj.P.Val <= pCutoff & abs(logFC) >= lfcCutoff)

str(sigDE)
# 'data.frame':	12915 obs. of  4 variables:

length(unique(sigDE$cds))
# [1] 6962

write.table(
  sigDE, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/dea_sigDE_stats.txt"
)

write(unique(sigDE$cds), "results/sigDE_cdsID.txt")
