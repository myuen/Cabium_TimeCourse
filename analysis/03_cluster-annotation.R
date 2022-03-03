library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)


### Add annotations to differentially expressed contigs from 3 sources:

### (1) Manually curated sequences from C. Keeling from interior Spruce PG-29
### (2) InterProScan.
### (3) BLAST RefSeq Plant


# ===


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


sigDE_stats <- left_join(sigDE_stats, clusters, by = "cds")

# Since we have contigs differentially expressed in multiple comparisons so we
# have a larger count of members in respective clusters.
table(sigDE_stats$cluster)
#    1    2    3    4    5 
# 2376 1737 2634 2443 3722 


### 0. Add clone ID that matches contigs
Clones <- 
  read.delim("data/clones.megablastAssembly.txt", 
             header = FALSE, stringsAsFactors = FALSE)

Clones <- Clones %>% select(V2, V1)

colnames(Clones) <- c("cds", "clone_id")



### 1. Process manual curate BLAST annotations

### Read BLAST results from PG29 manual annotation against assembly
Manual_Annot <- 
  read.delim("data/PG29ManualAnnot.blastpAssembly.txt",
             stringsAsFactors = FALSE, sep = "\t",
             header = FALSE)

# Header for tab-delimited BLAST results
colnames(Manual_Annot) <- 
  c("qseqid", "cds", "evalue", "bitscore",
    "length", "pident", "ppos", "qcovs", 
    "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "manual_annot")

str(Manual_Annot)
# 'data.frame':	2773 obs. of  16 variables:


# Only select query (i.e.contigs) with 95% or more identity and 80% to subject
# sequence (i.e. manual annotated sequences) length coverage and up-regulated
Manual_Annot <- Manual_Annot %>% 
  filter(cds %in% sigDE_stats$cds) %>% 
  filter(pident >= 95 & qcovs >= 80) %>% 
  select(cds, manual_annot)

str(Manual_Annot)
 # 'data.frame':	109 obs. of  2 variables:


#### 2. Process tsv InterProScan annotation
IPR <- read.delim("data/sigDE.aa.annotated.tsv.gz", 
                  stringsAsFactors = FALSE,
                  header = FALSE,
                  sep = "\t")

colnames(IPR) <- 
  c("cds", "MD5", "length", "analysis", 
    "sig_acc", "sig_desc",
    "start_loc", "stop_loc",
    "score",
    "status",
    "date",
    "ipr_annot", "ipr_annot_desc", "pathways")

str(IPR)
# 'data.frame':	56154 obs. of  14 variables:

# Cast score as a numeric
IPR$score <- as.numeric(IPR$score)


# Subset InterProScan annotation 
# Remove NA from coercion.
# Remove MobiDBLite annotations.  It does not add any value to our purpose.
# Retain only those with a e-value 1e-20 or below
IPR_Filtered <- IPR %>% 
  filter(!is.na(score)) %>% 
  filter(score <= 1e-20) %>% 
  filter(analysis != "MobiDBLite") %>% 
  select(cds, sig_acc, sig_desc, ipr_annot, ipr_annot_desc, sig_desc)


# Number of counts of contigs with InterProScan annotations
length(IPR_Filtered$cds)
# [1] 29924

# Number of counts of unique contigs with InterProScan annotations.  What is
# shown is there are a lot of contigs with multiple InterProScan match.
length(unique(IPR_Filtered$cds))
# [1] 5574


# Compress multi-line match
IPR_Filtered <- IPR_Filtered %>%
  group_by(cds) %>%
  nest()


IPR_Filtered$data <-
  map_df(IPR_Filtered$data, function(n){
    
    sig_accs <- unlist(n["sig_acc"])
    sig_accs <- sig_accs %>% 
      str_to_upper() %>% 
      unique()
    
    # Remove cells with no info
    sig_descs <- unlist(n["sig_desc"])
    sig_descs <- sig_descs[sig_descs != "-"] %>% 
      str_to_upper() %>% 
      unique()
  
    ipr_annots <- unlist(n["ipr_annot"])
    ipr_annots <- ipr_annots[ipr_annots != "-"] %>% 
      str_to_upper() %>% 
      unique()

    ipr_annot_descs <- unlist(n["ipr_annot_desc"])
    ipr_annot_descs <- ipr_annot_descs[ipr_annot_descs != "-"] %>% 
      str_to_upper() %>% 
      unique()
    
    # Collapse multi-line annotation into single line
    df <- data.frame(
      sig_accs <-str_c(sig_accs, collapse = ";"),
      
      sig_descs <- str_c(sig_descs, collapse = ";"),
      
      ipr_annots <- str_c(ipr_annots, collapse = ";"),

      ipr_annot_descs <- str_c(ipr_annot_descs, collapse = ";"))
})

IPR_Filtered <- IPR_Filtered %>% unnest(cols = c(data))

colnames(IPR_Filtered) <- c("cds", "ipr_sig_accs", "ipr_sig_descs",
                            "ipr_annots", "ipr_annot_descs")

str(IPR_Filtered)
# grouped_df [5,574 × 5] (S3: grouped_df/tbl_df/tbl/data.frame)



#### 3. Process BLAST tsv annotations
RefSeq_Annot <- 
  read.delim("data/sigDE.pep.blastpRefSeqPlant208.txt",
             stringsAsFactors = FALSE, sep = "\t",
             header = FALSE)

#### Header for tab-delimited BLAST results
colnames(RefSeq_Annot) <- 
  c("cds", "BLAST_hit", "evalue", "bitscore", "length", "pident", "ppos", 
    "qcovs", "qlen", "qstart", "qend", "qframe", 
    "slen", "sstart", "send", "sframe", 
    "staxid", "sskingdom", "sscinames", "scomnames", "BLAST_desc")

RefSeq_Annot <- RefSeq_Annot %>% 
  select(cds, BLAST_hit, BLAST_desc)

str(RefSeq_Annot)
# 'data.frame':	5773 obs. of  3 variables:


#### 4. Append all results into a final table
sigDE_Annot <- left_join(sigDE_stats, Clones, by = "cds")
sigDE_Annot <- left_join(sigDE_Annot, Manual_Annot, by = "cds")
sigDE_Annot <- left_join(sigDE_Annot, IPR_Filtered, by = "cds")
sigDE_Annot <- left_join(sigDE_Annot, RefSeq_Annot, by = "cds")

str(sigDE_Annot)
# 'data.frame':	13676 obs. of  12 variables:


write_delim(sigDE_Annot, "results/sigDE_annotations.txt.gz",
            delim = "\t", na = "")
