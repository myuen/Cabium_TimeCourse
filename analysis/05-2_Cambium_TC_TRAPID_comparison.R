library(dplyr)
library(stringr)

cambium_tc_structural <- 
  read.delim("data/TRAPID_comparison/cambium_tc_structural_data.txt")

str(cambium_tc_structural)
# 'data.frame':	47338 obs. of  11 variables:
  
colnames(cambium_tc_structural)[1] <- "transcript_id"


cambium_tc_protein_domain <- 
  read.delim("data/TRAPID_comparison/cambium_tc_protein_domain.txt")

cambium_tc_protein_domain <- cambium_tc_protein_domain |> select(-"X.counter")

str(cambium_tc_protein_domain)
# 'data.frame':	118953 obs. of  3 variables:

cambium_tc_fl <- cambium_tc_structural |> 
  filter(meta_annotation == "Full Length")

cambium_tc_fl <- left_join(cambium_tc_fl, cambium_tc_protein_domain,
                           by = "transcript_id")
str(cambium_tc_fl)
# 'data.frame':	60276 obs. of  13 variables:


### Transcript IDs
# BZIP
bzips <- cambium_tc_fl |> 
  filter(description == "Basic-leucine zipper domain")
str(bzips)
# 'data.frame':	44 obs. of  13 variables:

write.table(unique(bzips$transcript_id), "results/cambium_tc_bzip_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# HSF
hsfs <- cambium_tc_fl |> 
  filter(description == 'Heat shock transcription factor family')
str(hsfs)
# 'data.frame':	13 obs. of  13 variables:

write.table(unique(hsfs$transcript_id), "results/cambium_tc_hsfs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# MADS
mads <- cambium_tc_fl |> 
  filter(description == 'Transcription factor, MADS-box')
str(mads)
# 'data.frame':	38 obs. of  13 variables:

write.table(unique(mads$transcript_id), "results/cambium_tc_mads_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# NAC
nacs <- cambium_tc_fl |> 
  filter(description == 'NAC domain')
str(nacs)
# 'data.frame':	49 obs. of  13 variables:

write.table(unique(nacs$transcript_id), "results/cambium_tc_nacs_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

### Terpene Synthases
tps <- cambium_tc_fl |> 
  filter(str_detect(cambium_tc_fl$description, "Terpene synthase"))
str(tps)
# 'data.frame':	26 obs. of  13 variables:

write.table(unique(tps$transcript_id), "results/cambium_tc_tps_ids.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
