library(dplyr)
library(stringr)


### Of interest to us are two main categories:
### (1) gene involved in MEV/MEP pathway, terpene synthases and P450s.
### (2) Transcription factors that might be involved in the regulation of
### expression of genes identified in (1) as well as the development of
### traumatic resin duct in spruce cambium tissue.

#### Best sources for genes involved in (1) comes from PG-29 manual annotation
#### from C. Keeling. Transcription factors (TF) for (2) relies on annotation
#### from InterProScan

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

table(clusters$cluster)
#    1    2    3    4    5 
# 1234  998 1368 1396 1964 


sigDE <- left_join(sigDE_stats, clusters)

# Number of members in respective clusters
table(sigDE$cluster)
#    1    2    3    4    5 
# 2376 1737 2634 2443 3722 


#### Add BLAST annotations to result

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
sigDE_manualAnnot <- blastRslt %>% 
  filter(sseqid %in% sigDE$cds) %>% 
  filter(pident >= 95 & qcovs >= 80) %>% 
  select(sseqid, qseqid)

str(sigDE_manualAnnot)
# 'data.frame':	109 obs. of  2 variables:


# Add gene function to manual annotation
manualAnnotGeneFunc <- 
  read.delim("data/ManualAnnotationMatrix.txt",
             stringsAsFactors = FALSE, header = FALSE)

str(manualAnnotGeneFunc)
# 'data.frame':	565 obs. of  2 variables:

sigDE_manualAnnot <- 
  left_join(sigDE_manualAnnot, manualAnnotGeneFunc, 
            by = c("qseqid" = "V1"))

sigDE_manualAnnot <- sigDE_manualAnnot %>% 
  select(-qseqid) %>% 
  unique()

colnames(sigDE_manualAnnot) <- c("cds", "desc")

sigDE_manualAnnot$annot_method <- "manual"

str(sigDE_manualAnnot)
# 'data.frame':	99 obs. of  3 variables:



#### Add InterProScan annotation
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


# Subset InterProScan annotation
sigDE_IPR_Annot <- ipr %>% 
  filter(cds %in% sigDE$cds) %>%
  # filter(str_detect(tolower(ipr_annot_desc), "transcription factor")) %>% 
  select(cds, ipr_annot, ipr_annot_desc, sig_desc)

colnames(sigDE_IPR_Annot) <- c("cds", "annot_id", "desc")

#### WORK IN PROGRESS ### 
# There are multiple hits from InterProScan for 1 contigs.
# Append all results in a single line.

mDat <- sigDE_IPR_Annot %>% 
  filter(cds == "TRINITY_DN70060_c0_g1_i1.p1")

nDat <- mDat %>%
  group_by(cds) %>% 
  nest()

map_df(nDat$data, function(n){
  ipr_annots <- unlist(n["ipr_annot"])
  ipr_annots <- ipr_annots[ipr_annots != "-"]
  
  df <- data.frame(
    ipr_annots = str_c(ipr_annots, collapse = ";"),
    sig_desc = paste(unlist(n["sig_desc"]), collapse = ";"))
})


###

sigDE_IPR_Annot$annot_method <- "InterProScan"

all_annots <- rbind(sigDE_manualAnnot, sigDE_IPR_Annot)

str(all_annots)
# 'data.frame':	441 obs. of  4 variables:


sigDE_annot <- left_join(sigDE, all_annots)

str(sigDE_annot)
# 'data.frame':	13271 obs. of  8 variables:

sigDE_annot <- sigDE_annot[order(sigDE_annot$cluster),]

write.table(sigDE_up_annot, "results/sigDE_up_annotation.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")

