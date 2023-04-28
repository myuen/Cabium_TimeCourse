library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyr)

# Read standardized expression from Mfuzz
s.data <- read.table("results/standardized_expression.txt", 
                     header = TRUE, stringsAsFactors = FALSE)
str(s.data)
# 'data.frame':	6762 obs. of  6 variables:


# Read cluster members
clstr_members <- read.table("results/cluster_members.feb17.txt", 
                            header = TRUE, stringsAsFactors = FALSE)
str(clstr_members)
# 'data.frame':	6762 obs. of  2 variables:


# Transform table into long format and prepare for plot
s.data <- s.data %>% 
  pivot_longer(2:6, names_to = "library", values_to = "expression")

# Extract timepoint from library name
s.data$timepoint <- str_extract(s.data$library, "\\d$")

s.data$timepoint <- as.integer(s.data$timepoint)

# Convert timepoint to 'day'
s.data[s.data$timepoint == 1, 4] <- 0
s.data[s.data$timepoint == 2, 4] <- 1
s.data[s.data$timepoint == 3, 4] <- 2
s.data[s.data$timepoint == 4, 4] <- 4
s.data[s.data$timepoint == 5, 4] <- 8

 
# Join both data.frame
clstr_members_exprs <- left_join(s.data, clstr_members, by = 'cds')


# Read transcription factor contig IDs
bzip <- read.csv("data/bzip.ctgIDs.txt", stringsAsFactors = FALSE)
bzip <- bzip %>%
  add_column('enzyme' = 'tf_bzip')

hsf <- read.csv("data/hsf.ctgIDs.txt", stringsAsFactors = FALSE)
hsf <- hsf %>%
  add_column('enzyme' = 'tf_hsf')

mads <- read.csv("data/mads.ctgIDs.txt", stringsAsFactors = FALSE)
mads <- mads %>%
  add_column('enzyme' = 'tf_mads')

nac <- read.csv("data/nac.ctgIDs.txt", stringsAsFactors = FALSE)
nac <- nac %>%
  add_column('enzyme' = 'tf_nac')

tf <- bind_rows(bzip, hsf, mads, nac)
colnames(tf)[1] <- 'cds'

# Read MEP enzyme contig IDs
mep <- read.csv("data/MEP.txt", stringsAsFactors = FALSE, col.names = "cds")
mep <- mep %>% add_column('enzyme' = "mep")

# Read TPS enzyme contig IDs
tps <- read.csv("data/TPS.txt", stringsAsFactors = FALSE, col.names = "cds")
tps <- tps %>% add_column('enzyme' = "tps")


# Stack all tables
enzymes <- rbind(tf, mep, tps)

# Now the final table have expression, cluster number and enzyme type
clstr_members_exprs2 <- 
  left_join(clstr_members_exprs, enzymes, by = "cds")

# Replace NA as character string na
clstr_members_exprs2$enzyme <- 
  replace_na(clstr_members_exprs2$enzyme, "na")

# Create an 'alpha' column with default value as 1.  
# This column is for alpha value in plot
clstr_members_exprs2$alpha_value <- 1

# For cluster don't belong to any 50% of their expression value
clstr_members_exprs2[clstr_members_exprs2$enzyme == 'na', 'alpha_value'] <- 
  abs(clstr_members_exprs2[clstr_members_exprs2$enzyme == 'na', 'expression']) * 0.5

# Cast enzyme as factors
clstr_members_exprs2$enzyme <- 
  factor(clstr_members_exprs2$enzyme, 
         levels = c('tf_bzip', 'tf_hsf', 'tf_mads', 'tf_nac', 'mep', 'tps', 'na'))

clstr_members_exprs2 <- 
  clstr_members_exprs2[order(clstr_members_exprs2$enzyme),]


# Plot expression graph
g <-
  ggplot(clstr_members_exprs2,
         aes(x = timepoint, y = expression, group = rev(cds), 
             color = enzyme)) + 

  geom_line(aes(alpha = alpha_value)) + 
  
  # Break graph by cluster
  facet_wrap(. ~ cluster, nrow = 3, ncol = 2, 
              labeller = label_both) +
  
  # Colour code for different enzymes
  scale_color_manual(#na.value = 'grey85', 
    values = c('tf_bzip' = 'red',
               'tf_hsf' = 'red',
               'tf_mads' = 'red',
               'tf_nac' = 'red',
               'mep' = 'blue',
               'tps' = 'black',
               'na' = 'grey85')) +
  
  # X-axis break points
  scale_x_continuous(labels = c(0, 1, 2, 4, 8),
                     limits = c(0, 8),
                     breaks = c(0, 1, 2, 4, 8)) +
  
  # Plot, x-axis, y-axis label
  ggtitle("Soft Clustering of Differentially Expressed Contigs") +
  xlab("Time (days)") + 
  ylab("Normalized expression") +

  # Aesthetics 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey25'),
        legend.position = 'none', 
        axis.line = element_line(linetype = 1),
        strip.background = element_blank(), 
        strip.placement = "outside")

ggsave("results/figures/soft_clustering_plot.feb22.pdf", g, 
       # width = 2500, height = 2000, units = 'px',
       width = 8, height = 8, units = 'in',
       dpi = 320)

