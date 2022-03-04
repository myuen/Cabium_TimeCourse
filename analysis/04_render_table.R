library(dplyr)
library(gt)
library(readr)
library(stringr)
library(tibble)

annots <- read_tsv("results/sigDE_annotations.txt.gz", 
                     show_col_types = FALSE)

annots$logFC <- as.double(annots$logFC)

annots$adj.P.Val <- as.double(annots$adj.P.Val)

annots$ipr_sig_accs <- annots$ipr_sig_accs %>% 
  str_replace_all(";", ", ") %>% 
  str_wrap(width = 80)

annots$ipr_sig_descs <- 
  str_wrap(annots$ipr_sig_descs, width = 80)

annots$ipr_annots <- annots$ipr_annots %>% 
  str_replace_all(";", ", ") %>% 
  str_wrap(width = 80)

annots$ipr_annot_descs <- 
  str_wrap(annots$ipr_annot_descs, width = 80)

colnames(annots) <- c("CDS", "log2FC", "Adj. p-value", "Contrast", 
                      "Cluster", "Clone ID", "Manual Annot",
                      "IPR Sigature Acc", "IPR Signature Desc", 
                      "IPR Annot", "IPR Annot Desc",
                      "BLAST Top Hit", "BLAST Desc")

annots <- annots %>% 
  select("CDS", "Cluster", "Clone ID",
         "Contrast", "log2FC", "Adj. p-value",
         "Manual Annot",
         "IPR Sigature Acc", "IPR Signature Desc", 
         "IPR Annot", "IPR Annot Desc",
         "BLAST Top Hit", "BLAST Desc")

str(annots)
# spec_tbl_df [13,678 × 13] (S3: spec_tbl_df/tbl_df/tbl/data.frame)


### Table group by contigs
tbl_groupby_ctg <- annots %>%
  group_by(CDS) %>% 
  gt() %>% 

  # Header
  tab_header(title = md("**Annotations for All Differentially Expressed Contigs**")) %>%
  # Footer
  tab_source_note(source_note = md("**InterProScan and BLAST threshold at 1e-20**")) %>% 

  # Setting up subtitles  
  tab_spanner(
    label = md("**Differential Expression Stats**"),
    columns = c("log2FC", "Adj. p-value", "Contrast")) %>% 
  
  tab_spanner(
    label = md("**InterProScan**"),
    columns = c("IPR Sigature Acc", "IPR Signature Desc",
                "IPR Annot", "IPR Annot Desc")) %>% 

  tab_spanner(
    label = md("**BLAST RefSeq Plant**"),
    columns = c("BLAST Top Hit", "BLAST Desc")) %>% 

  cols_align(
    align = "center",
    columns = everything()) %>% 
  
  cols_align(
    align = "right",
    columns = "CDS")

gtsave(tbl_groupby_ctg, "results/annotaion.group_by_ctg.html")

#####

tbl_groupby_cluster <- annots %>%
  group_by(Cluster) %>%
  gt() %>%

  # Header
  tab_header(title = md("**Annotations for All Differentially Expressed Contigs**")) %>%
  # Footer
  tab_source_note(source_note = md("**InterProScan and BLAST threshold at 1e-20**")) %>%

  tab_spanner(
    label = md("**Differential Expression Stats**"),
    columns = c("log2FC", "Adj. p-value", "Contrast")) %>%

  tab_spanner(
    label = md("**InterProScan**"),
    columns = c("IPR Sigature Acc", "IPR Signature Desc",
                "IPR Annot", "IPR Annot Desc")) %>%

  tab_spanner(
    label = md("**BLAST RefSeq Plant**"),
    columns = c("BLAST Top Hit", "BLAST Desc")) %>% 
  
  cols_align(
    align = "center",
    columns = everything()) %>% 
  
  cols_align(
    align = "right",
    columns = "CDS")

gtsave(tbl_groupby_cluster, "results/annotaion.group_by_ctg.html")


### Group by Clone ID

tbl_groupby_clone <- annots %>%
  filter(`Clone ID` != "") %>% 
  select(-`Manual Annot`) %>% 
  group_by(`Clone ID`) %>%
  gt(rowname_col = "CDS") %>% 
  
  # Header
  tab_header(title = md("**Annotations for All Differentially Expressed Contigs**")) %>%
  # Footer
  tab_source_note(source_note = md("**InterProScan cutoff  at 1e-20**")) %>%
  tab_source_note(source_note = md("**BLAST threshold at 1e-20**")) %>%
  
  tab_spanner(
    label = md("**Differential Expression Stats**"),
    columns = c("log2FC", "Adj. p-value", "Contrast")) %>%
  
  tab_spanner(
    label = md("**InterProScan**"),
    columns = c("IPR Sigature Acc", "IPR Signature Desc",
                "IPR Annot", "IPR Annot Desc")) %>%
  
  tab_spanner(
    label = md("**BLAST RefSeq Plant**"),
    columns = c("BLAST Top Hit", "BLAST Desc")) %>% 

  cols_align(
    align = "center",
    columns = everything()) %>% 

  cols_align(
    align = "right",
    columns = "CDS")
  
# tbl_groupby_clone
gtsave(tbl_groupby_clone, "results/annotaion.group_by_clone.html")


### Group by transcription factors

tbl_groupby_tf <- annots %>% 
  filter(str_detect(`IPR Signature Desc`, "TRANSCRIPTION FACTOR")) %>% 
  gt() %>% 
  
  # Header
  tab_header(title = md("**Annotations for All Differentially Expressed Putative Transcript Factors**")) %>%
  # Footer
  tab_source_note(source_note = md("**InterProScan and BLAST threshold at 1e-20**")) %>%
  
  tab_spanner(
    label = md("**Differential Expression Stats**"),
    columns = c("log2FC", "Adj. p-value", "Contrast")) %>%
  
  tab_spanner(
    label = md("**InterProScan**"),
    columns = c("IPR Sigature Acc", "IPR Signature Desc",
                "IPR Annot", "IPR Annot Desc")) %>%
  
  tab_spanner(
    label = md("**BLAST RefSeq Plant**"),
    columns = c("BLAST Top Hit", "BLAST Desc")) %>% 
  
  cols_align(
    align = "center",
    columns = everything()) %>% 
  
  cols_align(
    align = "right",
    columns = "CDS")

gtsave(tbl_groupby_tf, "results/annotaion.group_by_tf.html")


### Group by manual annotation

tbl_groupby_ma <- annots[order(annots$log2FC, decreasing = TRUE),] %>% 
  filter(!is.na(`Manual Annot`)) %>% 
  group_by(`Manual Annot`) %>%
  gt(rowname_col = "CDS") %>%
  
  # Header
  tab_header(title = md("**Annotations for All Differentially Expressed Contigs**")) %>%
  # Footer
  tab_source_note(source_note = md("**InterProScan e-value threshold at 1e-20**")) %>%
  tab_source_note(source_note = md("**BLAST e-valuethreshold at 1e-20**")) %>%
  
  tab_spanner(
    label = md("**Differential Expression Stats**"),
    columns = c("log2FC", "Adj. p-value", "Contrast")) %>%
  
  tab_spanner(
    label = md("**InterProScan**"),
    columns = c("IPR Sigature Acc", "IPR Signature Desc",
                "IPR Annot", "IPR Annot Desc")) %>%
  
  tab_spanner(
    label = md("**BLAST RefSeq Plant**"),
    columns = c("BLAST Top Hit", "BLAST Desc")) %>% 
  
  cols_align(
    align = "right",
    columns = "CDS") %>% 
  
  cols_align(
    !starts_with("CDS"),
    align = "center",
    columns = everything()) %>% 
  
  cols_width(
    CDS ~ px(300),
    everything() ~ px(150)
  )


gtsave(tbl_groupby_ma, "results/annotaion.group_by_ma.html")

