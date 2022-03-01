library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

### Add annotations to differentially expressed contigs from 3 sources:

### (1) Manually curated sequences from C. Keeling from interior Spruce PG-29
### (2) InterProScan.
### (3) BLAST RefSeq Plant


#===


#### Read DE statistics
sigDE_stats <- read.delim("results/dea_sigDE_stats.txt",
                          stringsAsFactors = FALSE)

str(sigDE_stats)
# 'data.frame':	12915 obs. of  4 variables:



### Assign cluster ID to DE contigs
clusters <- read.table("results/cluster_members.txt",
                       header = TRUE)

str(clusters)
# 'data.frame':	6960 obs. of  2 variables:

# Number of members in respective clusters
table(clusters$cluster)
#    1    2    3    4    5 
# 1234  998 1368 1396 1964 


sigDE <- left_join(sigDE_stats, clusters)

# Since we have contigs differentially expressed in multiple comparisons so we
# have a larger count of members in respective clusters.
table(sigDE$cluster)
#    1    2    3    4    5 
# 2376 1737 2634 2443 3722 


#### Process manual curate BLAST annotations

# Read BLAST results from PG29 manual annotation against assembly
blastRslt <- 
  read.delim("data/PG29ManualAnnot.blastpAssembly.txt",
             stringsAsFactors = FALSE, sep = "\t",
             header = FALSE)

#### Header for tab-delimited BLAST results
colnames(blastRslt) <- 
  c("qseqid", "sseqid", "evalue", "bitscore",
    "length", "pident", "ppos", "qcovs", 
    "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "salltitles")

str(blastRslt)
# 'data.frame':	2773 obs. of  16 variables:


# Only select query (i.e.contigs) with 95% or more identity and 80% to subject
# sequence (i.e. manual annotated sequences) length coverage and up-regulated
sigDE_manual_Annot <- blastRslt %>% 
  filter(sseqid %in% sigDE$cds) %>% 
  filter(pident >= 95 & qcovs >= 80) %>% 
  select(sseqid, qseqid)

str(sigDE_manual_Annot)
 # 'data.frame':	109 obs. of  2 variables:


# Add gene function to manual annotation
manualAnnotGeneFunc <- 
  read.delim("data/ManualAnnotationMatrix.txt",
             stringsAsFactors = FALSE, header = FALSE)

str(manualAnnotGeneFunc)
# 'data.frame':	565 obs. of  2 variables:

sigDE_manual_Annot <- 
  left_join(sigDE_manual_Annot, manualAnnotGeneFunc, 
            by = c("qseqid" = "V1"))

sigDE_manual_Annot <- sigDE_manual_Annot %>% 
  select(-qseqid) %>% 
  unique()

colnames(sigDE_manual_Annot) <- c("cds", "manual annotation")

# sigDE_manual_Annot$annot_method <- "manual"

str(sigDE_manual_Annot)
# 'data.frame':	99 obs. of  3 variables:



#### Process tsv InterProScan annotation
ipr <- read.delim("data/sigDE.aa.annotated.tsv.gz", 
                  stringsAsFactors = FALSE,
                  header = FALSE,
                  sep = "\t")

colnames(ipr) <- 
  c("cds", "MD5", "length", "analysis", 
    "sig_acc", "sig_desc",
    "start_loc", "stop_loc",
    "score",
    "status",
    "date",
    "ipr_annot", "ipr_annot_desc", "pathways")

str(ipr)
# 'data.frame':	56154 obs. of  14 variables:

# Cast score as a numeric
ipr$score <- as.numeric(ipr$score)


# Subset InterProScan annotation 
# Remove NA from coercion.
# Retain only those with a e-value 1e-20 or below
ipr <- ipr %>% 
  filter(!is.na(score)) %>% 
  filter(score <= 1e-20) %>% 
  # filter()
  select(cds, ipr_annot, ipr_annot_desc, sig_desc)


# Number of counts of contigs with InterProScan annotations
length(sigDE_IPR_Annot$cds)
# [1] 56154

# Number of counts of unique contigs with InterProScan annotations.  What is
# shown is there are a lot of contigs with multiple InterProScan match.
length(unique(sigDE_IPR_Annot$cds))
# [1] 6537


# Compress multi-line match
sigDE_IPR_Annot <- sigDE_IPR_Annot %>%
  group_by(cds) %>% 
  nest()

sigDE_IPR_Annot$data <- 
  map_df(sigDE_IPR_Annot$data, function(n){

  annot_ids <- unlist(n["annot_id"])
  annot_ids <- annot_ids[annot_ids != "-"]

  ipr_annot_descs <- unlist(n["ipr_annot_desc"])
  ipr_annot_descs <- ipr_annot_descs[ipr_annot_descs != "-"]
  
  sig_descs <- unlist(n["sig_desc"])
  sig_descs <- sig_descs[sig_descs != "-"]

  
  df <- data.frame(
    annot_ids = str_c(annot_ids, collapse = ";"),
    ipr_annot_descs = str_c(ipr_annot_descs, collapse = ";"),
    sig_descs = str_c(sig_descs, collapse = ";"))
})

sigDE_IPR_Annot <- sigDE_IPR_Annot %>% unnest(cols = c(data))

sigDE_IPR_Annot$annot_method <- "InterProScan"

###


# all_annots <- rbind(sigDE_manual_Annot, sigDE_IPR_Annot)

# str(all_annots)
# 'data.frame':	441 obs. of  4 variables:


# Add manual annotation to DE contigs
sigDE_annot <- left_join(sigDE, sigDE_manual_Annot)
# Joining, by = "cds"

str(sigDE_annot)
# 'data.frame':	12956 obs. of  7 variables:


# For contigs without manual annotation, we will add InterProScan functional
# annotation

sigDE_annot <- sigDE_annot[order(sigDE_annot$cluster),]

write.table(sigDE_up_annot, "results/sigDE_up_annotation.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")

