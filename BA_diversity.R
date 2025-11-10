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

metadata1c<-readr::read_tsv("meema_baby_fecal_metadata.tsv")
View(metadata1b)

setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/Kiersten")
#Reading Taxonomy
taxonomy1c<-read_qza("MEEMA_B_Fec_taxonomy.qza")

#break up this string and for that purpose the parse_taxonomy() function is provided:
taxonomy1c<-parse_taxonomy(taxonomy1c$data)

#Creating a Phyloseq Object
MEEMA_B_Fec<-qza_to_phyloseq(
  features="MEEMA_B_Fec_table.qza",
  tree="MEEMA_B_Fec_rooted-tree.qza",
  taxonomy="MEEMA_B_Fec_taxonomy.qza",
  metadata = "meema_baby_fecal_metadata.tsv"
)

#View the phyloseq data 
MEEMA_B_Fec

MEEMA_B_Fec_filt = prune_taxa(phyloseq::taxa_sums(MEEMA_B_Fec) > 1, MEEMA_B_Fec)

# prune_taxa() alters your table based on specified parameters
# taxa_sums() refers to your bacterial species and their counts
# in my example I am removing samples do not appear at least twice in the table

#######BETA DIVERSITY ANALYSIS#############

# Create Distance Matrix & Data Frame
df_MEEMA_B <- as(sample_data(MEEMA_B_Fec_filt), "data.frame")

#Center log-ratio transformation - Step 1 of getting Aitchison distances
df_MEEMA_B_clr <- microbiome::transform(MEEMA_B_Fec_filt, "clr")

#Generate distance matrix (Aitchison) via Euclidean AFTER CLR
aitch_MEEMA_B <- phyloseq::distance(df_MEEMA_B_clr, method = "euclidean")

#weighted Unifrac
B_UniFrac_Weighted <- UniFrac(MEEMA_B_Fec_filt, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

#Regular PCoA for aitchison
MEEMA_B_PCA <- ordinate(MEEMA_B_Fec_filt, method = "PCoA", distance = "aitch_MEEMA_B") 

MEEMA_B_graph_aitch <- plot_ordination(MEEMA_B_Fec_filt, MEEMA_B_PCA, color = "Group", shape = "Timepoint") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse() +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  ggtitle("Aitchison Distance PCoA")

#PCoA for Weighted Unifrac
MEEMA_B_PCA_unifrac <- ordinate(MEEMA_B_Fec_filt, method = "PCoA", distance = "B_UniFrac_Weighted") 

MEEMA_B_graph2 <- plot_ordination(MEEMA_B_Fec_filt, MEEMA_B_PCA_unifrac, color = "Group", shape = "Timepoint") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse() +
  scale_color_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  ggtitle("Weighted UniFrac Distance PCoA")

#Redundancy analysis - Please fix this, axes are massive
MEEMA_B_PCA2 <- ordinate(MEEMA_B_Fec_filt, method = "RDA", distance = "aitch_MEEMA_B") 

MEEMA_B_graph3 <- plot_ordination(MEEMA_B_Fec_filt, MEEMA_B_PCA2, color = "Group") +
  theme_pubr() +
  labs(color = "", shape = "") +
  stat_ellipse()

#### Create Taxa Data Frame ######
MEEMA_B_Fec_filt
MEEMA_B_genera <- MEEMA_B_Fec_filt %>%
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

alpha.table.B1 <- microbiome::alpha(MEEMA_B_Fec_filt, index = "all")

# Create SampleID Column
alpha.table.B1$X.SampleID <- SummarizedExperiment::rownames(alpha.table.B1)
df_MEEMA_B$X.SampleID <- SummarizedExperiment::rownames(df_MEEMA_B)

# Merge Tables
alpha.table.B1 <- merge(alpha.table.B1, df_MEEMA_B, by = "X.SampleID")

# Plot
(ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_inverse_simpson, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Inverse Simpson index") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)

(pBM1b <- ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_inverse_simpson, fill = Group)) +
    geom_boxplot() +
    stat_smooth(method = "lm") +
    stat_cor() +
    labs(x = "Timepoint", y = "Inverse Simpson index") +
    scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    theme_minimal()
)

get_box_stats.B1 <- function(y, upper_limit = max(alpha.table.B1$diversity_inverse_simpson) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}


get_box_stats.B2.1 <- function(y, upper_limit = max(alpha.table.B1$diversity_shannon) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}


get_box_stats.B3.1 <- function(y, upper_limit = max(alpha.table.B1$observed) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

(ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_shannon, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Shannon Diversity index") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)

(ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_shannon, fill = Group)) +
    geom_boxplot() +
    stat_smooth(method = "lm") +
    stat_cor() +
    labs(x = "Timepoint", y = "Shannon Diversity index") +
    scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    theme_minimal()
)


ggplot(alpha.table.B1, aes(x = Timepoint, y = observed, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1,60))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats.B2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = -5.2) +
  labs(x = "Timepoint", y = "Observed ASVs") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()

ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_inverse_simpson, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,12))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats.B2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = -3.5) +
  labs(x = "Timepoint", y = "Inverse Simpson Index") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()

ggplot(alpha.table.B1, aes(x = Timepoint, y = diversity_shannon, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1.5,3))+
  stat_compare_means(fontface = "bold") +
  stat_summary(fun.data = get_box_stats.B2.1, geom = "text", 
               position = position_dodge(width = 0.75),
               hjust = 0.5, vjust = 0.75) +
  labs(x = "Timepoint", y = "Shannon Index Diversity") +
  scale_fill_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
  theme_minimal()


(ggplot(alpha.table.B1, aes(x = Timepoint, y = observed, color = Group, group = Baseline.cpartid)) +
    geom_point(size = 3.5) + 
    geom_line(aes(linetype=Group), size = 1.2, show.legend = TRUE) +
    scale_linetype_manual(values = c("Control" = "solid", "Intervention" = "dashed")) +
    labs(x = "Timepoint", y = "Observed ASVs") +
    scale_colour_manual(values = c("Control" = "#9c3b5e", "Intervention" = "#e191af")) +
    scale_x_discrete(expand = c(0.09, 0)) + scale_y_continuous(expand = c(0.15, 0)) +
    theme_minimal()
)
