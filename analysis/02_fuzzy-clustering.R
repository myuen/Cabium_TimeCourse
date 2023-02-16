library(Mfuzz)

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyr)
library(tidyselect)


# Read normalized CPM
cpm <- read.table("results/tmm_normalized.cpm.txt")


# Function to calculate row variance
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


cpm.stats <-
  cpm %>% transmute(
    "MJ.T1" = rowMeans(select(cpm, starts_with("MJ.T1"))),
    "MJ.T2" = rowMeans(select(cpm, starts_with("MJ.T2"))),
    "MJ.T3" = rowMeans(select(cpm, starts_with("MJ.T3"))),
    "MJ.T4" = rowMeans(select(cpm, starts_with("MJ.T4"))),
    "MJ.T5" = rowMeans(select(cpm, starts_with("MJ.T5"))),
    
    "TW.T1" = rowMeans(select(cpm, starts_with("TW.T1"))),
    "TW.T2" = rowMeans(select(cpm, starts_with("TW.T2"))),
    "TW.T3" = rowMeans(select(cpm, starts_with("TW.T3"))),
    "TW.T4" = rowMeans(select(cpm, starts_with("TW.T4"))),
    "TW.T5" = rowMeans(select(cpm, starts_with("TW.T5"))),
    
    "MJ.mean" = rowMeans(select(cpm, starts_with("M"))),
    "TW.mean" = rowMeans(select(cpm, starts_with("T"))),
    
    "MJ.var" = rowVar(select(cpm, starts_with("M"))),
    "TW.var" = rowVar(select(cpm, starts_with("T")))
  )


# Read statistics for differential expression sequences
sigDE <- read.table("results/dea_sigDE_stats.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE)

str(sigDE)
# 'data.frame':	12915 obs. of  4 variables:


# Filter CPM if it is sig. differentially expressed
sigDE.stats <- 
  cpm.stats[rownames(cpm.stats) %in% sigDE$cds,] %>% 
  filter(MJ.var != 0) %>% 
  filter(TW.var != 0) %>% 
  select(1:5)

str(sigDE.stats)
# 'data.frame':	6762 obs. of  5 variables:


# Add time to table
time <- c(2, 24, 48, 96, 192)

sigDE.stats <- rbind(time, sigDE.stats)

rownames(sigDE.stats)[1] <- "time"


# table2eset function must read from a file
write.table(sigDE.stats, file = "results/sigDE.stats.tbl.new.txt",
            sep = '\t', quote = F, col.names = NA)


# Convert an expression table to eset object
data <- table2eset("results/sigDE.stats.tbl.new.txt")


# Standardise data
data.s <- standardise(data)


(m1 <- mestimate(data.s))
# [1] 2.004545

# pdf("results/figures/dmin.pdf", width = 20, height = 10)

Dmin(data.s, m = m1, crange = seq(2, 10, 1), repeats = 5, visu = TRUE)
# [1] 2.2942950 2.3378309 2.0339764 1.1626455 1.1139902 1.0279514 0.8745872 0.7924840 0.6685793

dev.off()


clust <- 5

c <- mfuzz(data.s, c = clust, m = m1)

# pdf("results/figures/dea-clusters.new.pdf", width = 20, height = 10)
# 
# mfuzz.plot2(data.s, cl = c, mfrow = c(2, 3), colo = "fancy",
#             time.labels = c(2, 24, 48, 96, 192),
#             xlab = "Time (hrs)",
#             ylab = expression("Log"[2]*'FC'),
#             centre = TRUE,
#             # Increase font size of title
#             cex.main = 1.2,
#             # Increase font size of tick point
#             cex.axis = 1.1,
#             # Increase font size of axis lable
#             cex.lab = 1.2,
#             x11 = FALSE)
# 
# dev.off()


#####

### Standardize expression data
s.data <- exprs(data.s) %>% 
  as.data.frame() %>% 
  rownames_to_column('cds')


# Transform table into long format and prepare for plot
s.data <- s.data %>% 
  pivot_longer(2:6, names_to = "library", values_to = "expression")

# extract timepoint from library name
s.data$timepoint <- str_extract(s.data$library, "\\d$")

s.data$timepoint <- as.integer(s.data$timepoint)

# Convert timepoint to 'hours'
s.data[s.data$timepoint == 5, 4] <- 192
s.data[s.data$timepoint == 4, 4] <- 96
s.data[s.data$timepoint == 3, 4] <- 48
s.data[s.data$timepoint == 2, 4] <- 24
s.data[s.data$timepoint == 1, 4] <- 2


# Write out cluster members
clstr_members <- as.data.frame(c[3])

clstr_members <- clstr_members %>% 
  rownames_to_column("cds")

# write.table(clstr_members, "results/cluster_members.txt", 
#             quote = FALSE, sep = "\t", row.names = FALSE)

write.table(clstr_members, "results/cluster_members.new.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

##### 

clstr_members_exprs <- left_join(s.data, clstr_members, by = 'cds')

#####

bzip <- read.csv("data/TRAPID_TF_BZIP_FL.csv", stringsAsFactors = FALSE)
bzip <- bzip %>%
  select(2) %>%
  add_column('tf' = 'bzip')

hsf <- read.csv("data/TRAPID_TF_HSF_FL.csv", stringsAsFactors = FALSE)
hsf <- hsf %>%
  select(2) %>%
  add_column('tf' = 'hsf')

mads <- read.csv("data/TRAPID_TF_MADS_FL.csv", stringsAsFactors = FALSE)
mads <- mads %>%
  select(2) %>%
  add_column('tf' = 'mads')

nac <- read.csv("data/TRAPID_TF_NAC_FL.csv", stringsAsFactors = FALSE)
nac <- nac %>%
  select(2) %>%
  add_column('tf' = 'nac')


tf <- bind_rows(bzip, hsf, mads, nac)
colnames(tf)[1] <- 'cds'

clstr_members_exprs2 <- left_join(clstr_members_exprs, tf, by = "cds")

clstr_members_exprs2$tf <- factor(clstr_members_exprs2$tf, 
                                  levels = c('bzip', 'hsf', 'mads', 'nac'))

#####

jpeg("results/figures/draft.jpg", width = 3000, height = 1600, units = 'px')
# pdf("results/figures/draft.pdf", width = 20, height = 10)

# clstr_members_exprs2$tf <- replace_na(clstr_members_exprs2$tf, "na")
clstr_members_exprs2 <- clstr_members_exprs2[order(clstr_members_exprs2$tf),]

# clstr_members_exprs2$tf <- replace_na(clstr_members_exprs2$tf, "na")


g <-
  ggplot(clstr_members_exprs2,
       aes(x = timepoint, y = expression, group = rev(cds), color = tf)) + 
  geom_line() +
  facet_wrap( ~ cluster) + 
  
  scale_color_manual(na.value = 'grey90', values = c('bzip' = 'blue',
                                'hsf' = 'coral',
                                'mads' = 'olivedrab',
                                'nac' = 'red')) +
    
  # scale_linewidth_manual(values = c('bzip' = 1,
  #                                   'hsf' = 1,
  #                                   'mads' = 1,
  #                                   'nac' = 1,
  #                                   'na' = 0.5)) +

  scale_alpha_manual(na.value = 0.1, values = c('bzip' = 1,
                                'hsf' = 1,
                                'mads' = 1,
                                'nac' = 1)) + 

  scale_x_continuous(labels = c(2, 24, 48, 96, 192),
                     limits = c(2, 192),
                     breaks = c(2, 24, 48, 96, 192)) +
  
  xlab("Time Point (Hours)") + 
  ylab("Expression") +
  ggtitle("Soft Clustering of differential expression contigs") +
  
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

ggsave("results/figures/draft.jpg", g, width = 3000, height = 1500, units = 'px')

# last_plot() + aes(group=rev(cds))

dev.off()

