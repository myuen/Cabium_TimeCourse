library(dplyr)
library(ggplot2)
library(tidyr)


plot.data <- data.frame()


### BZIP -----------------------------------------------------------------------
crd.bzip <- scan("results/crd_bzip_ids.txt", what='character')
trd.bzip <- scan("results/cambium_tc_bzip_ids.txt", what='character')

bzip.crd_vs_trd <- 
  read.table("data/TRAPID_comparison/bzip.crd_blastp_cambium.txt")
bzip.trd_vs_crd <- 
  read.table("data/TRAPID_comparison/bzip.cambium_blastp_crd.txt")

# V6 is percent identity and V8 percent coverage
bzip.crd_vs_trd <- 
  bzip.crd_vs_trd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V1, V2)
colnames(bzip.crd_vs_trd) <- c('crd_id', 'trd_id')
str(bzip.crd_vs_trd)
# 'data.frame':	35 obs. of  2 variables:

bzip.trd_vs_crd <- 
  bzip.trd_vs_crd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V2, V1)
colnames(bzip.trd_vs_crd) <- c('crd_id', 'trd_id')
str(bzip.trd_vs_crd)
# 'data.frame':	33 obs. of  2 variables:

bzip.common <-
  dim(intersect(bzip.crd_vs_trd, bzip.trd_vs_crd))[1]

bzip <- data.frame(
  'tf' = 'bzip',
  'crd' = length(crd.bzip) - bzip.common,
  'common' = bzip.common,
  'trd' = length(trd.bzip) - bzip.common)

plot.data <- rbind(plot.data, bzip)


### HSFS -----------------------------------------------------------------------
crd.hsfs <- scan("results/crd_hsfs_ids.txt", what='character')
trd.hsfs <- scan("results/cambium_tc_hsfs_ids.txt", what='character')

hsfs.crd_vs_trd <- 
  read.table("data/TRAPID_comparison/hsfs.crd_blastp_cambium.txt")

# V6 is percent identity and V8 percent coverage
hsfs.crd_vs_trd <- 
  hsfs.crd_vs_trd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V1, V2)
colnames(hsfs.crd_vs_trd) <- c('crd_id', 'trd_id')
str(hsfs.crd_vs_trd)
# 'data.frame':	6 obs. of  2 variables:


hsfs.trd_vs_crd <- 
  read.table("data/TRAPID_comparison/hsfs.cambium_blastp_crd.txt")

hsfs.trd_vs_crd <- 
  hsfs.trd_vs_crd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V2, V1)
colnames(hsfs.trd_vs_crd) <- c('crd_id', 'trd_id')
str(hsfs.trd_vs_crd)
# 'data.frame':	8 obs. of  2 variables:

hsfs.common <-
  dim(intersect(hsfs.crd_vs_trd, hsfs.trd_vs_crd))[1]

hsfs <- data.frame(
  'tf' = 'hsfs',
  'crd' = length(crd.hsfs) - hsfs.common,
  'common' = bzip.common,
  'trd' = length(trd.hsfs) - hsfs.common)

plot.data <- rbind(plot.data, hsfs)


### MADS -----------------------------------------------------------------------

crd.mads <- scan("results/crd_mads_ids.txt", what='character')
trd.mads <- scan("results/cambium_tc_mads_ids.txt", what='character')

mads.crd_vs_trd <- 
  read.table("data/TRAPID_comparison/mads.crd_blastp_cambium.txt")
mads.trd_vs_crd <- 
  read.table("data/TRAPID_comparison/mads.cambium_blastp_crd.txt")

# V6 is percent identity and V8 percent coverage
mads.crd_vs_trd <- 
  mads.crd_vs_trd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V1, V2)
colnames(mads.crd_vs_trd) <- c('crd_id', 'trd_id')
str(mads.crd_vs_trd)
# 'data.frame':	35 obs. of  2 variables:

mads.trd_vs_crd <- 
  mads.trd_vs_crd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V2, V1)
colnames(mads.trd_vs_crd) <- c('crd_id', 'trd_id')
str(mads.trd_vs_crd)
# 'data.frame':	33 obs. of  2 variables:

mads.common <-
  dim(intersect(mads.crd_vs_trd, mads.trd_vs_crd))[1]

mads <- data.frame(
  'tf' = 'mads',
  'crd' = length(crd.mads) - mads.common,
  'common' = mads.common,
  'trd' = length(trd.mads) - mads.common)

plot.data <- rbind(plot.data, mads)


### NACS -----------------------------------------------------------------------
crd.nac <- scan("results/crd_nacs_ids.txt", what='character')
trd.nac <- scan("results/cambium_tc_nacs_ids.txt", what='character')

nac.crd_vs_trd <- 
  read.table("data/TRAPID_comparison/nacs.crd_blastp_cambium.txt")
nac.trd_vs_crd <- 
  read.table("data/TRAPID_comparison/nacs.cambium_blastp_crd.txt")

# V6 is percent identity and V8 percent coverage
nac.crd_vs_trd <- 
  nac.crd_vs_trd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V1, V2)
colnames(nac.crd_vs_trd) <- c('crd_id', 'trd_id')
str(nac.crd_vs_trd)
# 'data.frame':	17 obs. of  2 variables:

nac.trd_vs_crd <- 
  nac.trd_vs_crd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V2, V1)
colnames(nac.trd_vs_crd) <- c('crd_id', 'trd_id')
str(nac.trd_vs_crd)
# 'data.frame':	14 obs. of  2 variables:

nac.common <-
  dim(intersect(nac.crd_vs_trd, nac.trd_vs_crd))[1]

nac <- data.frame(
  'tf' = 'nacs',
  'crd' = length(crd.nac) - nac.common,
  'common' = nac.common,
  'trd' = length(trd.nac) - nac.common)

plot.data <- rbind(plot.data, nac)


### TPS ------------------------------------------------------------------------
crd.tps <- scan("results/crd_tps_ids.txt", what='character')
trd.tps <- scan("results/cambium_tc_tps_ids.txt", what='character')

tps.crd_vs_trd <- 
  read.table("data/TRAPID_comparison/tps.crd_blastp_cambium.txt")
tps.trd_vs_crd <- 
  read.table("data/TRAPID_comparison/tps.cambium_blastp_crd.txt")

# V6 is percent identity and V8 percent coverage
tps.crd_vs_trd <- 
  tps.crd_vs_trd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V1, V2)
colnames(tps.crd_vs_trd) <- c('crd_id', 'trd_id')
str(tps.crd_vs_trd)
# 'data.frame':	2 obs. of  2 variables:

tps.trd_vs_crd <- 
  tps.trd_vs_crd |> 
  filter(V6 >= 99 & V8 >= 99) |> 
  select(V2, V1)
colnames(tps.trd_vs_crd) <- c('crd_id', 'trd_id')
str(tps.trd_vs_crd)
# 'data.frame':	2 obs. of  2 variables:

tps.common <-
  dim(intersect(tps.crd_vs_trd, tps.trd_vs_crd))[1]

tps <- data.frame(
  'tf' = 'tps',
  'crd' = length(crd.tps) - tps.common,
  'common' = tps.common,
  'trd' = length(trd.tps) - tps.common)

plot.data <- rbind(plot.data, tps)


### House keeping genes --------------------------------------------------------

crd.house_keeping <- 
  read.table("data/TRAPID_comparison/house-keeping-genes.blastp_crd.txt")
trd.house_keeping <- 
  read.table("data/TRAPID_comparison/house-keeping-genes.blastp_trd.txt")

# V6 is percent identity and V8 percent coverage
crd.house_keeping <-
  crd.house_keeping |>
  filter(V6 >= 90 & V8 >= 90) |>
  select(V1, V2)
colnames(crd.house_keeping) <- c('ath_id', 'crd_id')
str(crd.house_keeping)
# 'data.frame':	4 obs. of  2 variables:

trd.house_keeping <-
  trd.house_keeping |>
  filter(V6 >= 90 & V8 >= 90) |>
  select(V1, V2)
colnames(trd.house_keeping) <- c('ath_id', 'trd_id')
str(trd.house_keeping)
# 'data.frame':	4 obs. of  2 variables:

# Find house keeping genes that are expressed in both crd and trd
house_keeping.common <-
  length(intersect(crd.house_keeping$ath_id, trd.house_keeping$ath_id))


house_keeping <- data.frame(
  'tf' = 'house keeping',
  'crd' = dim(crd.house_keeping)[1] - house_keeping.common,
  'common' = house_keeping.common,
  'trd' = dim(trd.house_keeping)[1] - house_keeping.common)

plot.data <- rbind(plot.data, house_keeping)


### Plot -----------------------------------------------------------------------

plot.data.long <- 
  plot.data |> 
  pivot_longer(cols = 2:4, names_to = "type", values_to = "count")

plot.data.long$tf <- toupper(plot.data.long$tf)

# Replace zero counts to NA for text label in plot
plot.data.long$count <- na_if(plot.data.long$count, 0)

plot.data.long$type <-
  factor(plot.data.long$type, levels=c('trd', 'common', 'crd'))

# Set the max for x-axis
x_max <- plot.data.long |> 
  group_by(tf) |> 
  summarise(sum=sum(count)) |> 
  select(sum) |> 
  max()

x_max <- round(x_max, digits = -1)


(g <- ggplot(plot.data.long, aes(y=tf, x=count)) + 
    geom_col(width=0.5, aes(fill= factor(type))) +

    # Add label
    geom_text(aes(label = count), 
              colour = 'white',
              family = 'helvetica',
              fontface = 'bold',
              position = position_stack(vjust=0.5)) +
    
    # Plot title and axis label
    labs(title = "Number of transcripton factors and terpene synthase expressed") +
    xlab("Count") + 
    ylab("") +
    
    # Add tick mark on x-axis
    scale_x_continuous(breaks=seq(0, x_max, 10)) +
    scale_y_discrete(limits=c("HOUSE KEEPING", "TPS", "NACS", "MADS", "HSFS", "BZIP")) +


    # Manual edit colors
    scale_fill_manual(breaks=c('crd', 'common', 'trd'), 
                      labels = c('crd' = 'CRD-only', 
                                 'common' = 'Both', 
                                 'trd' = 'TRD-only'),
                      values = c('#748d9a', '#CE9916', '#AD000B')) +
    
    theme_minimal() +
  
    # Remove legend title
    theme(legend.title = element_blank())
)
