library(dplyr)
library(stringr)

crd_structural <- 
  read.delim("data/TRAPID_comparison/crd_structural_data.txt")

str(crd_structural)
# 'data.frame':	47445 obs. of  11 variables:

colnames(crd_structural)[1] <- "transcript_id"


crd_protein_domain <- 
  read.delim("data/TRAPID_comparison/crd_protein_domain.txt")

crd_protein_domain <- crd_protein_domain |> select(-"X.counter")

str(crd_protein_domain)
# 'data.frame':	119732 obs. of  3 variables:


crd_fl <- crd_structural |> filter(meta_annotation == "Full Length")

crd_fl <- left_join(crd_fl, crd_protein_domain, 
                     by = "transcript_id")

str(crd_fl)
# 'data.frame':	59921 obs. of  13 variables:


### Transcript IDs
# BZIP
bzips <- crd_fl |> filter(description == "Basic-leucine zipper domain")
str(bzips)
# 'data.frame':	47 obs. of  13 variables:

write.table(unique(bzips$transcript_id), "results/crd_bzip_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# HSF
hsfs <- crd_fl |> filter(description == 'Heat shock transcription factor family')
str(hsfs)
# 'data.frame':	13 obs. of  13 variables:

write.table(unique(hsfs$transcript_id), "results/crd_hsfs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# MADS
mads <- crd_fl |> 
  filter(description == 'Transcription factor, MADS-box')
str(mads)
# 'data.frame':	30 obs. of  13 variables:

write.table(unique(mads$transcript_id), "results/crd_mads_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# NAC
nacs <- crd_fl |> 
  filter(description == 'NAC domain')
str(nacs)
# 'data.frame':	33 obs. of  13 variables:

write.table(unique(nacs$transcript_id), "results/crd_nacs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

### Terpene Synthases
tps <- 
  crd_fl |> filter(str_detect(crd_fl$description, "Terpene synthase"))
str(tps)
# 'data.frame':	68 obs. of  13 variables:

write.table(unique(tps$transcript_id), "results/crd_tps_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
