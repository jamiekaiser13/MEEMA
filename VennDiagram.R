#Making a Venn Diagram for the overlap between taxonomy for the moms'#
#breastmilk and the babys' feces#

library(knitr)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(ggforce)
library(grid)
library(qiime2R)
library(dplyr)
library(stringr)
library(ggvenn)
library(ggVennDiagram)
install.packages("VennDiagram")
library(VennDiagram)

#create table with taxonomy
View(inverted_table_2024)

#making dictionary
taxonomy_jamie_2024_dict <- setNames(taxonomy_jamie_2024$data$Taxon, taxonomy_jamie_2024$data$Feature.ID) 
#removing species from dictionary
taxonomy_jamie_2024_togenus <- str_remove(taxonomy_jamie_2024$data$Taxon, "; s__.*")
taxonomy_jamie_2024_dict_togenus <- setNames(taxonomy_jamie_2024_togenus, taxonomy_jamie_2024$data$Feature.ID)

#rename feature tables with taxonomy
colnames(table_2024) <- taxonomy_jamie_2024_dict_togenus[colnames(inverted_table_2024)]
View(table_2024)

ma_rows <- grepl("^MA", rownames(table_2024))
ma_rows
mb_rows <- grepl("^MB", rownames(table_2024))
mb_rows
ba_rows <- grepl("^BA", rownames(table_2024))
ba_rows

###remove MA samples
features_ba_mb_2024 <- table_2024[!ma_rows, ]
View(features_ba_mb_2024)
#remove columns (taxa) that are not present in MB or BA samples
features_ba_mb_2024_filt <- features_ba_mb_2024[,colSums(features_ba_mb_2024) > 0]
View(features_ba_mb_2024_filt)

###BA samples
features_ba_2024 <- features_ba_mb_2024_filt[ba_rows, ]
View(features_ba_2024)
#remove columns (taxa) that are not present in BA samples
features_ba_2024_filt <- features_ba_2024[,colSums(features_ba_2024) > 0]
View(features_ba_2024_filt)
#making a row that sums the totals for the BA samples
features_ba_2024_sums <- features_ba_2024_filt %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum)))
View(features_ba_2024_sums)
baby_sums <- features_ba_2024_sums["...14", ]
rownames(baby_sums) <- "BA_total"
View(baby_sums)
baby_taxa <- colnames(baby_sums)

#MB_samples
mb_rows_filt <- grepl("^MB", rownames(features_ba_mb_2024_filt))
features_mb_2024 <- features_ba_mb_2024_filt[mb_rows_filt, ]
View(features_mb_2024)
#remove columns (taxa) that are not present in MB  samples
features_mb_2024_filt <- features_mb_2024[,colSums(features_mb_2024) > 0]
View(features_mb_2024_filt)
#making a row that sums the totals for the MB samples
features_mb_2024_sums <- features_mb_2024_filt %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum)))
View(features_mb_2024_sums)
mb_sums <- features_mb_2024_sums["...14", ]
rownames(mb_sums) <- "MB_total"
View(mb_sums)
mb_taxa <- colnames(mb_sums)

#making venn diagram
taxa_mb_ba <- list(baby_venn = baby_taxa,
                   mb_venn = mb_taxa)
mb_ba_venn <- ggVennDiagram(taxa_mb_ba,
                            label = "percent", 
                            label_size = 8,
                            category.names = c("baby taxa", "breastmilk taxa")) +
  scale_fill_gradient(low = "#d3b5d3", high = "#786178") +
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "Taxa overlap of baby gut and mother breastmilk microbiomes")
print(mb_ba_venn)

#getting the values in the overlap
mb_ba_partitions <- get.venn.partitions(taxa_mb_ba)
View(mb_ba_partitions)
#cleaning the data up so we can actually evaluate the overlap
colnames(mb_ba_partitions) <- c("baby_venn", "mb_venn", "set", "values", "count")
mb_ba_overlap <- mb_ba_partitions[1, "values"]
mb_ba_overlap
#retaining only genus to help identify taxa
extract_genus_regex <- function(taxonomy_string) {
  result <- sub(".*?(g__.*)$", "\\1", taxonomy_string)
  ifelse(grepl("g__", taxonomy_string), result, NA)
}
mb_ba_overlap_genus <- sapply(mb_ba_overlap, extract_genus_regex)
mb_ba_overlap_genus


#combining MB and BA sums to use if needed
mb_ba_sums <- rbind(baby_sums, mb_sums)
View(mb_ba_sums)
mb_ba_sums[, 1:5]


