library(dplyr)
library(stringr)

bark_structural <- 
  read.delim("data/TRAPID_comparison/bark_lmd_structural_data.txt")

str(bark_structural)
# 'data.frame':	58189 obs. of  11 variables:

colnames(bark_structural)[1] <- "transcript_id"


bark_protein_domain <- 
  read.delim("data/TRAPID_comparison/bark_lmd_protein_domain.txt")

bark_protein_domain <- bark_protein_domain |> select(-"X.counter")

str(bark_protein_domain)
# 'data.frame':	142386 obs. of  3 variables:


bark_fl <- bark_structural |> filter(meta_annotation == "Full Length")

bark_fl <- left_join(bark_fl, bark_protein_domain, 
                     by = "transcript_id")

str(bark_fl)
# 'data.frame':	58419 obs. of  13 variables:

### Transcript IDs
# BZIP
bzips <- bark_fl |> filter(description == "Basic-leucine zipper domain")
str(bzips)
# 'data.frame':	41 obs. of  13 variables:

write.table(unique(bzips$transcript_id), "results/bark_lmd_bzip_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# HSF
hsfs <- bark_fl |> filter(description == 'Heat shock transcription factor family')
str(hsfs)
# 'data.frame':	14 obs. of  13 variables:

write.table(unique(hsfs$transcript_id), "results/bark_lmd_hsfs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# MADS
mads <- bark_fl |> 
  filter(description == 'Transcription factor, MADS-box')
str(mads)
# 'data.frame':	28 obs. of  13 variables:

write.table(unique(mads$transcript_id), "results/bark_lmd_mads_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# NAC
nacs <- bark_fl |> 
  filter(description == 'NAC domain')
str(nacs)
# 'data.frame':	41 obs. of  13 variables:

write.table(unique(nacs$transcript_id), "results/bark_lmd_nacs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

### Terpene Synthases
tps <- 
  bark_fl |> filter(str_detect(bark_fl$description, "Terpene synthase"))
str(tps)
# 'data.frame':	38 obs. of  13 variables:

write.table(unique(tps$transcript_id), "results/bark_lmd_tps_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
