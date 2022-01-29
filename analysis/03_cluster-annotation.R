library(dplyr)
library(stringr)

### Of interest to us are two main categories:
### (1) gene involved in MEV/MEP pathway, terpene synthases and P450s.
### (2) Transcription factors that might be involved in the regulation of
### expression of genes identified in (1) as well as the development of
### traumatic resin duct in spruce cambium tissue.

#### Best source for genes involved in (1) comes from PG-29 manual annotation
#### from C. Keeling. Transcription factors (TF) relied on annotation from
#### InterProScan

#===

#### Read DE statistics
sigDE_stats <- read.delim("results/dea_sigDE_stats.txt",
                          stringsAsFactors = FALSE)

str(sigDE_stats)
# 'data.frame':	12915 obs. of  4 variables:


#### Isolate sequence that are up-regulated
sigDE_up <- sigDE_stats %>% filter(logFC >= 0)

str(sigDE_up)
# 'data.frame':	7031 obs. of  4 variables:


### Assign cluster ID to DE contigs
clusters <- read.table("results/cluster_members.txt",
                       header = TRUE)

str(clusters)
# 'data.frame':	6960 obs. of  2 variables:

table(clusters$cluster)
#    1    2    3    4    5 
# 1396 1368 1234 1963  999 

sigDE_up <- left_join(sigDE_up, clusters)

# Number of members in respective clusters that are up-regulated
table(sigDE_up$cluster)
#    1    2    3    4    5 
# 1947  314 2310 1906  554 


#### Add BLAST annotations to result

#### Header for tab-delimited BLAST results

# Read BLAST results from PG29 manual annotation against assembly
blastRslt <- 
  read.delim("data/PG29ManualAnnot.blastpAssembly.txt",
             stringsAsFactors = FALSE, sep = "\t",
             header = FALSE)

colnames(blastRslt) <- 
  c("qseqid", "sseqid", "evalue", "bitscore",
    "length", "pident", "ppos", "qcovs", 
    "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "salltitles")

str(blastRslt)
# 'data.frame':	2773 obs. of  16 variables:

# Only select query (i.e.contigs) with 95% or more identity and 80% to subject
# sequence (i.e. manual annotated sequences) length coverage and up-regulated
sigDE_up_manualAnnot <- blastRslt %>% 
  filter(sseqid %in% sigDE_up$cds) %>% 
  filter(pident >= 95 & qcovs >= 80) %>% 
  select(sseqid, qseqid)

str(sigDE_up_manualAnnot)
# 'data.frame':	78 obs. of  2 variables:


# Add gene function to manual annotation
manualAnnotGeneFunc <- 
  read.delim("data/ManualAnnotationMatrix.txt",
             stringsAsFactors = FALSE, header = FALSE)

str(manualAnnotGeneFunc)
# 'data.frame':	565 obs. of  2 variables:

sigDE_up_manualAnnot <- 
  left_join(sigDE_up_manualAnnot, manualAnnotGeneFunc, 
            by = c("qseqid" = "V1"))

colnames(sigDE_up_manualAnnot) <- c("cds", "annot_id", "desc")

sigDE_up_manualAnnot$annot_method <- "manual"


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


# Subset InterProScan annotation, limit only to up-regulated contigs and only
# to TFs by keyword 'transcription factor'
sigDE_up_autoAnnot <- ipr %>% 
  filter(cds %in% sigDE_up$cds) %>%
  filter(str_detect(tolower(ipr_annot_desc), "transcription factor")) %>% 
  select(cds, ipr_annot_desc, sig_desc)

colnames(sigDE_up_autoAnnot) <- c("cds", "annot_id", "desc")

sigDE_up_autoAnnot$annot_method <- "InterProScan"

all_annots <- rbind(sigDE_up_manualAnnot, sigDE_up_autoAnnot)

str(all_annots)
# 'data.frame':	287 obs. of  4 variables:

sigDE_up_annot <- left_join(sigDE_up, all_annots)

str(sigDE_up_annot)
# 'data.frame':	7288 obs. of  8 variables:

sigDE_up_annot <- sigDE_up_annot[order(sigDE_up_annot$cluster),]

write.table(sigDE_up_annot, "results/sigDE_up_annotation.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
