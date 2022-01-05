library(dplyr)
library(stringr)

clusters <- read.table("results/cluster_members.txt", 
                       header = TRUE)

str(clusters)
# 'data.frame':	6960 obs. of  2 variables:

table(clusters$cluster)
#    1    2    3    4    5 
# 1963 1234 1368 1396  999 


ipr <- read.delim("data/sigDE.aa.annotated.tsv", 
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


ipr <- ipr %>% select(cds, length, analysis, 
                      sig_acc, sig_desc,
                      ipr_annot, ipr_annot_desc, pathways)

clusters.annot <- left_join(clusters, ipr)

clusters.annot$sig_desc <- str_to_upper(clusters.annot$sig_desc)

# Count number of transcription factor by cluster
clusters.annot %>% 
  filter(str_detect(sig_desc, "TRANSCRIPTION")) %>%
  group_by(cluster) %>% 
  select(cds, cluster) %>%
  unique() %>%
  count(cluster)

#   cluster     n
# 1       1   115
# 2       2    75
# 3       3    77
# 4       4    63
# 5       5    69

