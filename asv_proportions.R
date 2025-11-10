library(phyloseq)
install.packages("DESeq2")

#Importing feature table and taxonomy for the samples
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq")
features_jamie_2024 <- read_qza("table.qza")
taxonomy_jamie_2024 <- read_qza("taxonomy_jamie_2024.qza")
features_jamie_2024_df <- features_jamie_2024$data

#making dictionary
taxonomy_jamie_2024_dict <- setNames(taxonomy_jamie_2024$data$Taxon, taxonomy_jamie_2024$data$Feature.ID)
rm(features_jamie_2024)

#rename feature tables with taxonomy
row.names(features_jamie_2024_df) <- taxonomy_jamie_2024_dict[row.names(features_jamie_2024_df)]
View(features_jamie_2024_df)
write.csv(features_jamie_2024_df, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq/feature_table.csv")

feature_table_ba <- feature_table_baby$data
View(feature_table_ba)
row.names(feature_table_ba) <- taxonomy_jamie_2024_dict[row.names(feature_table_ba)]
View(feature_table_ba)
write.csv(feature_table_ba, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq/feature_table_baby.csv")

install.packages("readr")
library(readr)

bab_taxa_totals <- readr::read_tsv("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq/baby_taxa_totals.txt")
baby_taxa_totals <- bab_taxa_totals

install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
library(tidyverse)

baby_taxa_totals_positions <- baby_taxa_totals %>%
  mutate(csum = rev(cumsum(rev(`Total ASVs`))),
         pos = `Total ASVs`/2 + lead(csum, 1),
         pos = if_else(is.na(pos), `Total ASVs`/2, pos))

ggplot(baby_taxa_totals, aes(x='', y=`Total ASVs`, fill= fct_inorder(Genus))) + 
         geom_col() +
         coord_polar(theta = "y") +
         geom_label_repel(data = baby_taxa_totals_positions,
                   aes(y=pos, label = paste0(Genus)),
                   size = 4, nudge_x = 1, show.legend = FALSE) +
         guides(fill = guide_legend(title = "Genus")) + theme_void()

baby_taxa_sorted <- baby_taxa_totals[order(-baby_taxa_totals$`Total ASVs`),]
View(baby_taxa_sorted)  
baby_taxa_top10 <- baby_taxa_sorted[1:10,]
View(baby_taxa_top10)

baby_taxa_top10_positions <- baby_taxa_top10 %>%
  mutate(csum = rev(cumsum(rev(`Total ASVs`))),
         pos = `Total ASVs`/2 + lead(csum, 1),
         pos = if_else(is.na(pos), `Total ASVs`/2, pos))

ggplot(baby_taxa_top10, aes(x='', y=`Total ASVs`, fill= fct_inorder(Genus))) + 
  geom_col() +
  coord_polar(theta = "y") +
  geom_label_repel(data = baby_taxa_top10_positions,
                   aes(y=pos, label = paste0(Genus)),
                   size = 4, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Genus")) + theme_void()

baby_id <- c("BA001B", "BA001E", "BA002B", "BA002E", "BA003B", "BA003E", "BA005B", "BA005E", "BA006B", "BA007B", "BA008B")
features_ba <- features_jamie_2024$data[,baby_id]
View(features_ba)
features_ba <- data.frame(features_ba[rowSums(features_ba[])>0,])
View(features_ba)
typeof(features_ba)
  #[1] "list"
features_ba <- tibble::rownames_to_column(features_ba, "ASV")
View(features_ba)
features_ba$taxa <- taxonomy_jamie_2024_dict[features_ba$ASV]
View(features_ba)

feature_table_baby_csv <- read.csv("feature_table_baby.csv", stringsAsFactors = TRUE)
View(feature_table_baby_csv)
features_BA001B <- feature_table_baby_csv[,1:2]
View(features_BA001B)

features_BA001B_top10 <- features_BA001B[1:10, ]

#for loop to generate feature tables for each baby
ba_features_list <- list()
for (x in 2:9) {
  df_subname <- colnames(features_ba[x])
  df_name <- paste0("features_", df_subname)
  ba_features_list <- append(ba_features_list, df_name)
  assign(df_name, features_ba[ , c(1,x,13)])
}
print(ba_features_list)

ba_intervention_ids <- c("BA001B", "BA001E", "BA002B", "BA002E", "BA003B", "BA003E", "BA005B", "BA005E")

#rename the column that has the ASV counts to be called "count" in all dfs
for (x in ba_features_list) {
    df_name <- x
    df <- get(df_name)
    names(df)[2] <- "count"
    assign(x, df)
}
View(features_BA001B)

#remove rows that have a value of 0 in the count column
for(x in ba_features_list) {
  df_name <- x
  df <- get(df_name)
  df <- filter(df, `count` > 0)
  assign(x, df)
}
View(features_BA001B)

extract_genus_regex <- function(taxonomy_string) {
  result <- sub(".*?(g__.*)$", "\\1", taxonomy_string)
  ifelse(grepl("g__", taxonomy_string), result, NA)
}

for (x in ba_features_list) {
  df_name <- x
  df <- get(df_name)
  df2_name <- paste0(df_name, "_genus")
  df2 <- df %>% 
    mutate(genus = extract_genus_regex(taxa))
  assign(df2_name, df2)
}

for (x in ba_features_list) {
  df_name <- paste0(x, "_genus")
  print(df_name)
  df <- get(df_name)
  df2 <- df %>%
    mutate(genus = str_replace(df$genus, "g__; s__", NA_character_))
  assign(df_name, df2)
}
View(features_BA001B_genus)

for (x in ba_features_list) {
  df_name <- paste0(x, "_genus")
  df <- get(df_name)
  df$log_count <- log(df$count)
  assign(df_name, df)
}

for (x in ba_features_list) {
  df_name <- paste0(x, "_genus")
  df <- get(df_name)
  sum_count <- sum(df$count)
  df$relative_count <- (df$count / sum_count) * 100
  assign(df_name, df)
}

for (x in ba_features_list) {
  df_name <- paste0(x, "_genus")
  df <- get(df_name)
  df$genus <- sub(";.*", "", df$genus)
  assign(df_name, df)
}

features_BA001B_genus$log_count <- log(features_BA001B_genus$count)

ggplot(data = features_BA001B_genus, aes(x=genus, y=relative_count, fill = ASV)) + 
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) + 
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("BA001B ASVs")

library(stringr)
sum(str_detect(features_BA001B_genus$genus, "g__; s__"))


#trying to graph all together
features_BA001B_genus_labeled <- features_BA001B_genus  %>% 
  mutate(timepoint = "Baseline", sample_id = "BA001B")
features_BA001E_genus_labeled <- features_BA001E_genus  %>% 
  mutate(timepoint = "Endpoint", sample_id = "BA001E")
features_BA002B_genus_labeled <- features_BA002B_genus  %>% 
  mutate(timepoint = "Baseline", sample_id = "BA002B")
features_BA002E_genus_labeled <- features_BA002E_genus  %>% 
  mutate(timepoint = "Endpoint", sample_id = "BA002E")
features_BA003B_genus_labeled <- features_BA003B_genus  %>% 
  mutate(timepoint = "Baseline", sample_id = "BA003B")
features_BA003E_genus_labeled <- features_BA003E_genus  %>% 
  mutate(timepoint = "Endpoint", sample_id = "BA003E")
features_BA005B_genus_labeled <- features_BA005B_genus  %>% 
  mutate(timepoint = "Baseline", sample_id = "BA005B")
features_BA005E_genus_labeled <- features_BA005E_genus  %>% 
  mutate(timepoint = "Endpoint", sample_id = "BA005E")

genuses <- bind_rows(features_BA001B_genus_labeled, features_BA001E_genus_labeled,
                     features_BA002B_genus_labeled, features_BA002E_genus_labeled,
                     features_BA003B_genus_labeled, features_BA003E_genus_labeled,
                     features_BA005B_genus_labeled, features_BA005E_genus_labeled)

plot_by_genus <- ggplot(genuses, aes(x = timepoint, y = relative_count)) + 
  geom_boxplot(aes(fill = timepoint)) +
  facet_wrap(~genus, scales = "free_y") + 
  labs(y = "relative abundance", 
       title = "Relative Taxon Abundance in Baby Fecal Samples of Intervention Group") +
  theme(axis.title.x = element_blank(), 
      plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c("Baseline" = "#ccddff", "Endpoint" = "#6699ff"))

plot_by_genus

#graphing each baby baseline vs. endpoint

baby1 <- bind_rows(features_BA001B_genus_labeled, features_BA001E_genus_labeled)
baby2 <- bind_rows(features_BA002B_genus_labeled, features_BA002E_genus_labeled)
baby3 <- bind_rows(features_BA003B_genus_labeled, features_BA003E_genus_labeled)
baby5 <- bind_rows(features_BA005B_genus_labeled, features_BA005E_genus_labeled)

genuses_baby1 <- ggplot(baby1, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8, angle = -45)) +
  ggtitle("Relative Taxon Abundance in Baby 1 Fecal Samples")
genuses_baby1

genuses_baby2 <- ggplot(baby2, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8, angle= 90)) +
  ggtitle("Relative Taxon Abundance in Baby 2 Fecal Samples")
genuses_baby2

genuses_baby3 <- ggplot(baby3, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8)) +
  ggtitle("Relative Taxon Abundance in Baby 3 Fecal Samples")
genuses_baby3

genuses_baby5 <- ggplot(baby5, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8, angle = -45)) +
  ggtitle("Relative Taxon Abundance in Baby 5 Fecal Samples")
genuses_baby5


#percentages
genuses_baby1_percent <- ggplot(baby1, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "fill", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8)) +
  ggtitle("Relative Taxon Abundance in Baby 1 Fecal Samples")
genuses_baby1_percent

genuses_baby2_percent <- ggplot(baby2, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "fill", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8, angle= 90)) +
  ggtitle("Relative Taxon Abundance in Baby 2 Fecal Samples")
genuses_baby2_percent

genuses_baby3_percent <- ggplot(baby3, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "fill", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8)) +
  ggtitle("Relative Taxon Abundance in Baby 3 Fecal Samples")
genuses_baby3_percent

genuses_baby5_percent <- ggplot(baby5, aes(x = timepoint, y = relative_count, fill = ASV)) +
  geom_bar(position = "fill", stat = "identity", show.legend = FALSE) +
  facet_grid(~genus, switch = "x") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), 
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"),
        strip.text = element_text(size = 8, angle = -45)) +
  ggtitle("Relative Taxon Abundance in Baby 5 Fecal Samples")
genuses_baby5_percent


ggplot(genuses, aes(x = ASV, y = relative_count, fill = timepoint)) + 
  geom_boxplot() +
  facet_wrap(~genus, scales = "free") + 
  labs(y = "relative abundance", 
       title = "Relative Taxon Abundance in Baby Fecal Samples of Intervention Group") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c("Baseline" = "#ccddff", "Endpoint" = "#6699ff"))
