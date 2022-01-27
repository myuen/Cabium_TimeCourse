library(dplyr)
library(stringr)

# Reading and merging the information from different sources.


# Read DE statistics
sigDE_stats <- read.delim("results/dea_sigDE_stats.txt",
                          stringsAsFactors = FALSE)

str(sigDE_stats)
# 'data.frame':	12915 obs. of  4 variables:

sigDE_up <- sigDE_stats %>% filter(logFC >= 0)

str(sigDE_up)
# 'data.frame':	7031 obs. of  4 variables:


### Assign sigDE contigs to clusters
clusters <- read.table("results/cluster_members.txt",
                       header = TRUE)

str(clusters)
# 'data.frame':	6960 obs. of  2 variables:

table(clusters$cluster)
#    1    2    3    4    5 
# 1396 1368 1234 1963  999 

sigDE_up <- left_join(sigDE_up, clusters)


# Add BLAST results to annotation

# Header for tab-delimited BLAST results
blast_header <- c("qseqid", "sseqid", "evalue", "bitscore",
                  "length", "pident", "ppos", "qcovs", 
                  "qlen", "qstart", "qend", "qframe",
                  "slen", "sstart", "send", "salltitles")


# Read BLAST results from PG29 manual annotation against assembly
annot2de <- read.delim("data/PG29ManualAnnot.blastpAssembly.txt",
                       stringsAsFactors = FALSE, sep = "\t",
                       header = FALSE)

str(annot2de)
# 'data.frame':	2773 obs. of  16 variables:

colnames(annot2de) <- blast_header

# Only select query with 95% or more identity and 80% sequence length coverage
annot2de <- annot2de %>% filter(pident >= 95 & qcovs >= 80)

str(annot2de)
# 'data.frame':	238 obs. of  16 variables:


# Read BLAST results from assembly BLAST against PG29 manual annotations
de2annot <- read.delim("data/sigDE.blastpPG29ManualAnnot.txt", 
                     stringsAsFactors = FALSE, sep = "\t",
                     header = FALSE)

colnames(de2annot) <- blast_header

str(de2annot)
# 'data.frame':	2064 obs. of  16 variables:


# Find reciprocal best hit from BLAST results
reciprocalBest <- 
  inner_join(de2annot, annot2de, by = c("qseqid" = "sseqid", "sseqid" = "qseqid"))

str(reciprocalBest)
# 'data.frame':	104 obs. of  30 variables:

reciprocalBest <- reciprocalBest %>% 
  select(qseqid, sseqid)


# Function of manual annotated genes
geneFunc <- read.delim("data/ManualAnnotationMatrix.txt",
                       stringsAsFactors = FALSE, header = FALSE)

str(geneFunc)
# 'data.frame':	565 obs. of  2 variables:

sigDE_annotated <- 
  left_join(reciprocalBest, geneFunc, by = c("sseqid" = "V1"))

colnames(sigDE_annotated) <- c("cds", "gene", "function")

str(sigDE_annotated)
# 'data.frame':	104 obs. of  3 variables:

# write.table(sigDE_annotated, "results/sigDE_manual_annotation.txt",
            # col.names = TRUE, row.names = FALSE, 
            # sep = "\t", quote = FALSE)




upReg <- left_join(upReg, clusters)

str(upReg)
# 'data.frame':	7031 obs. of  5 variables:


ipr <- read.delim("data/sigDE.aa.annotated.tsv.gz", 
                  stringsAsFactors = FALSE,
                  header = FALSE,
                  sep = "\t")

str(ipr)
# 'data.frame':	56154 obs. of  14 variables:

colnames(ipr) <- 
  c("cds", "MD5", "length", "analysis", 
    "sig_acc", "sig_desc",
    "start_loc", "stop_loc",
    "score",
    "status",
    "date",
    "ipr_annot", "ipr_annot_desc", "pathways")


# Only keep selected columns
ipr <- ipr %>% 
  select(cds, length, analysis, 
         sig_acc, sig_desc,
         ipr_annot, ipr_annot_desc, pathways)


clusters.annot <- left_join(upReg, ipr)

clusters.annot$sig_desc <- str_to_upper(clusters.annot$sig_desc)

# Count number of transcription factor by cluster
clusters.annot %>% 
  filter(str_detect(sig_desc, "TRANSCRIPTION")) %>%
  group_by(cluster) %>% 
  select(cds, cluster) %>%
  unique() %>%
  count(cluster)

#   cluster     n
#     <int> <int>
# 1       1    47
# 2       2    72
# 3       3    60
# 4       4    66
# 5       5    37

# Write out transcription factor CDS IDs
tfs <- clusters.annot %>% 
  filter(str_detect(sig_desc, "TRANSCRIPTION")) %>%
  select(cds) %>%
  unique()

write(tfs$cds, "results/TF.cdsID.txt")
