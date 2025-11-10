#######################################################################
################## Mother's Fecal Analysis ############################
#######################################################################

outdated_packages_list <- list(c("Biobase", "Biostrings", "cli", "data.table",
                                 "IRanges", "multtest", "Rcpp", "rhdf5", "rhdf5filters",
                                 "Rhdf5lib", "rlang", "S4Vectors", "tibble",
                                 "vegan", "XVector", "zlibbioc"))
for(pkg in outdated_packages_list) {
  BiocManager::install(pkg)
}


###LOAD LIBRARIES###
library(rbiom)
library("gplots")
library(tidyverse)
library(qiime2R)
library("ggVennDiagram")
library(ggvenn)
library(ggplot2)
library(vegan)
library(package = "lattice")
library(permute)
library(dplyr)
library(magrittr)
library(scales)
library(ShortRead)
library(phyloseq)
library(qiime2R)
library(conflicted)
library(MetaTopics)
library(devtools)
library(usethis)
library(biomformat)
library(phyloseq)
library(ggpubr)
library(mixOmics)
library(ANCOMBC)
library(DT)
library(dplyr)
library(kableExtra)
library(nlme)
library(vegan)
library(DECIPHER)
library(phangorn)
library(gghighlight)
library(R.utils)

library(Biostrings)

library(microbiome)
library(ggsci)
library(reshape2)
library(phyloseq)
library(rbiom)
library(randomForestSRC)

library(RColorBrewer)

suppressPackageStartupMessages(library(microViz))

#Reading Metadata
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/Qiime2")

metadata1<-readr::read_tsv("meema_mom_fecal_metadata.tsv")
View(metadata1)

setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/Kiersten")
#Reading Taxonomy
taxonomy1<-read_qza("MEEMA_M_Fec_taxonomy.qza")

#break up this string and for that purpose the parse_taxonomy() function is provided:
taxonomy1<-parse_taxonomy(taxonomy1$data)

#Creating a Phyloseq Object
MEEMA_M_Fec<-qza_to_phyloseq(
  features="MEEMA_M_Fec_table.qza",
  tree="MEEMA_M_Fec_rooted-tree.qza",
  taxonomy="MEEMA_M_Fec_taxonomy.qza",
  metadata = "meema_mom_fecal_metadata.tsv"
)

#View the phyloseq data 
MEEMA_M_Fec

MEEMA_M_Fec_filt = prune_taxa(phyloseq::taxa_sums(MEEMA_M_Fec) > 1, MEEMA_M_Fec)

# prune_taxa() alters your table based on specified parameters
# taxa_sums() refers to your bacterial species and their counts
# in my example I am removing samples do not appear at least twice in the table

#######BETA DIVERSITY ANALYSIS#############

# Create Distance Matrix & Data Frame
df_MEEMA_M_Fec <- as(sample_data(MEEMA_M_Fec_filt), "data.frame")

#Center log-ratio transformation - Step 1 of getting Aitchison distances
d_MEEMA_M_Fec_clr <- microbiome::transform(MEEMA_M_Fec_filt, "clr")

#Generate distance matrix (Aitchison) via Euclidean AFTER CLR
aitch_MEEMA_M_Fec <- phyloseq::distance(d_MEEMA_M_Fec_clr, method = "euclidean")

#weighted Unifrac
M_Fec_UniFrac_Weighted <- UniFrac(MEEMA_M_Fec, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

#Regular PCoA for aitchison
MEEMA_M_Fec_PCA <- ordinate(MEEMA_M_Fec_filt, method = "PCoA", distance = "aitch_MEEMA_M_Fec") 

MEEMA_M_Fec_graph_aitch <- plot_ordination(MEEMA_M_Fec_filt, MEEMA_M_Fec_PCA, color = "Group", shape = "Timepoint") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse() +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))
  #ggtitle("Aitchison Distance PCoA")
MEEMA_M_Fec_graph_aitch

#PCoA for Weighted Unifrac
MEEMA_M_Fec_PCA_unifrac <- ordinate(MEEMA_M_Fec_filt, method = "PCoA", distance = "M_Fec_UniFrac_Weighted") 

MEEMA_M_Fec_graph2 <- plot_ordination(MEEMA_M_Fec_filt, MEEMA_M_Fec_PCA_unifrac, color = "Group", shape = "Timepoint") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse() +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  ggtitle("Weighted UniFrac Distance PCoA")
MEEMA_M_Fec_graph

#Redundancy analysis - Please fix this, axes are massive
MEEMA_M_Fec_PCA2 <- ordinate(MEEMA_M_Fec_filt, method = "RDA", distance = "aitch_MEEMA_M_Fec") 
MEEMA_M_Fec_PCA2

MEEMA_M_Fec_graph2 <- plot_ordination(MEEMA_M_Fec_filt, MEEMA_M_Fec_PCA2, color = "Group") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse()

MEEMA_M_Fec_graph2

#### Create Taxa Data Frame ######
MEEMA_M_Fec_filt
MEEMA_M_Fec_genera <- MEEMA_M_Fec_filt %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  dplyr::arrange(Phylum)                                      # Sort data frame alphabetically by phylum

##########################
### Diversity Analysis ###
##########################

### Alpha Diversity Analysis ###

# Alpha diversity measurees the number of distinct species within a a microbial
# group. There are an array of measures that look at very specific diversity
# indicies:
# - Shannon diversity: # of evenly distributed taxa
# - Simpson diversity: # of dominant taxa
# - Chao1: # of novel organisms in the community
# - Observed: the overall number of unique organisms in community
#
# Each has its own justiciation for measuring and analysis and is dependent
# on the question you are asking regarding your microbiome. I like to use the R
# package "microbiome" to assess alpha diversity. 
# 
# https://microbiome.github.io/tutorials/

alpha.table1 <- microbiome::alpha(MEEMA_M_Fec_filt, index = "all")

# Create SampleID Column
alpha.table1$X.SampleID <- SummarizedExperiment::rownames(alpha.table1)
df_MEEMA_M_Fec$X.SampleID <- SummarizedExperiment::rownames(df_MEEMA_M_Fec)

# Merge Tables
alpha.table1 <- merge(alpha.table1, df_MEEMA_M_Fec, by = "X.SampleID")
View(alpha.table1)

# Plot
(p1 <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_inverse_simpson, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Inverse Simpson index") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)

(p1b <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_inverse_simpson, fill = Group)) +
    geom_boxplot() +
    stat_smooth(method = "lm") +
    stat_cor() +
    labs(x = "Timepoint", y = "Inverse Simpson index") +
    scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    theme_minimal()
)

get_box_stats.1 <- function(y, upper_limit = max(alpha.table1$diversity_inverse_simpson) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}


get_box_stats2.1 <- function(y, upper_limit = max(alpha.table1$diversity_shannon) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}


get_box_stats3.1 <- function(y, upper_limit = max(alpha.table1$observed) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

(p2 <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_shannon, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Shannon Diversity index") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)

(p2b <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_shannon, fill = Group)) +
    geom_boxplot() +
    stat_smooth(method = "lm") +
    stat_cor() +
    labs(x = "Timepoint", y = "Shannon Diversity index") +
    scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    theme_minimal()
)


p3 <- ggplot(alpha.table1, aes(x = Timepoint, y = observed, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1,500))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = -4.5) +
  labs(x = "Timepoint", y = "Observed ASVs") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()

p4 <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_inverse_simpson, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1,80))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = -4.5) +
  labs(x = "Timepoint", y = "Inverse Simpson Index") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()

p6 <- ggplot(alpha.table1, aes(x = Timepoint, y = diversity_shannon, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(3,5.25))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = 0.3) +
  labs(x = "Timepoint", y = "Shannon Index Diversity") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()


ggplot(alpha.table1, aes(x = Group, y = observed, fill = Group)) +
  geom_point() +
  geom_boxplot(notch = TRUE, notchwidth = 0.93) +
  stat_compare_means() +
  stat_summary(fun.data = get_box_stats3, geom = "text", hjust = 0.5, vjust = 0.9) +
  labs(x = "Group type", y = "Observed ASVs") +
  theme_minimal()

(p5 <- ggplot(alpha.table1, aes(x = Timepoint, y = observed, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Observed ASVs") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)


#Diversity Metrics Against VO2Max
setwd("C:/Users/jmktw/OneDrive - University of Illinois Chicago/MEEMA Paper")

#uploading csv with MA metrics
deltas <- read.csv("Delta_Metrics.csv")

ggplot(deltas, aes(x=Delta_VO2Max, y=Delta_ASVs, color = Group)) +
  geom_point(size = 3.5) +
  labs(color = "") +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))

ggplot(deltas, aes(x=Delta_VO2Max, y=Delta_Shannon, color = Group)) +
  geom_point(size = 3.5) +
  labs(color = "") +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))

ggplot(deltas, aes(x=Delta_VO2Max, y=Delta_Simpson, color = Group)) +
  geom_point(size = 3.5) +
  labs(color = "") +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))

ggplot(deltas, aes(x=Delta_VO2Max, y=Aitchison, color = Group)) +
  geom_point(size = 3.5) +
  labs(color = "") +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))

ggplot(deltas, aes(x=Delta_VO2Max, y=Weighted_Unifrac, color = Group)) +
  geom_point(size = 3.5) +
  labs(color = "") +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af"))
