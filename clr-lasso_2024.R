#doing clr-lasso with the taxonomy generated from the 2024 version of qiime2

library(knitr)
library(glmnet)
library(selbal)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(ggforce)
library(grid)
library(qiime2R)
library(dplyr)
library(stringr)

#metadata file is the same as 2020 workflow, but I will have to upload a new feature table and taxonomy

#Importing feature table and taxonomy for the samples
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq")
features_jamie_2024 <- read_qza("table.qza")
taxonomy_jamie_2024 <- read_qza("taxonomy_jamie_2024.qza")

#inverting table
inverted_table_2024 <- data.frame(t(features_jamie_2024$data), check.names = FALSE)
View(inverted_table_2024)

#making dictionary
taxonomy_jamie_2024_dict <- setNames(taxonomy_jamie_2024$data$Taxon, taxonomy_jamie_2024$data$Feature.ID) 
#removing species from dictionary
taxonomy_jamie_2024_togenus <- str_remove(taxonomy_jamie_2024$data$Taxon, "; s__.*")
taxonomy_jamie_2024_dict_togenus <- setNames(taxonomy_jamie_2024_togenus, taxonomy_jamie_2024$data$Feature.ID)

#Prepping the feature data by adding 1 and taking the log of all features
features_jamie_2024_1 <- inverted_table_2024 + 1
features_jamie_2024_log <- log(features_jamie_2024_1)
View(features_jamie_2024_log)

#rename feature tables with taxonomy
colnames(features_jamie_2024_1) <- taxonomy_jamie_2024_dict_togenus[colnames(features_jamie_2024_1)]
colnames(features_jamie_2024_log) <- taxonomy_jamie_2024_dict_togenus[colnames(features_jamie_2024_log)]
View(features_jamie_2024_log)

##CLR-LASSO##
##Baby fecal samples of intervention group##

#restricting log feature table to only be the intervention baby samples
baby_intervention_id <- c("BA001B", "BA001E", "BA002B", "BA002E", "BA003B", "BA003E", "BA005B", "BA005E")
table_intervention_log_ba_2024 <- features_jamie_2024_log[baby_intervention_id, ]
View(table_intervention_log_ba_2024)
table_intervention_1_ba_2024 <- features_jamie_2024_1[baby_intervention_id, ]

#Clr function
clr_intervention_ba_2024 <- apply(table_intervention_log_ba_2024, 2, function(x) x- rowMeans(table_intervention_log_ba_2024))

#Determining the penalization parameter (lambda)
intervention_ba_2024.test_clrlasso <- glmnet(x = clr_intervention_ba_2024, y = ba_intervention, family = "binomial", nlambda = 50)
  # Warning message:
  #   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
  #               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_ba_2024.test_clrlasso, xvar = "lambda", label = T)
intervention_ba_2024.test_clrlasso
#A lambda of 0.02546 explains approximately 94% of variation with 5 features

intervention_ba_2024_clrlasso <- glmnet(x = clr_intervention_ba_2024, y = ba_intervention, family = "binomial", lambda = 0.03152)
  # Warning message:
  #   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
  #               one multinomial or binomial class has fewer than 8  observations; dangerous ground
intervention_ba_2024.results_clrlasso <- glmnet_wrapper(intervention_ba_2024_clrlasso, X = clr_intervention_ba_2024)
intervention_ba_2024.results_clrlasso$numVarSelect
#when I set the lambda to 0.02546, I get 9 features, which is not what the test output said it would be. Idk why. 
#When I set the lambda to 0.03152 (same lambda as when I used the 2020 data), I get 6 features, 
#which is not what the test output said it would be, but seems sufficient
#and similar to when I ran this function with the 2020 feature table. 
intervention_ba_2024.results_clrlasso$varSelect
# [1] "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__"                  
# [2] "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__europaeus"
# [3] "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__uniformis"        
# [4] "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__[Exiguobacteraceae]; g__; s__"                                 
# [5] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Dialister; s__"                      
# [6] "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales; f__; g__; s__" 

#selbal like plot for clr-lasso
intervention_ba_2024.clr_pos <- intervention_ba_2024.results_clrlasso$posCoefSelect
intervention_ba_2024.clr_nef <- intervention_ba_2024.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ba_2024.clr_pos), 
                 neg.names = names(intervention_ba_2024.clr_nef),
                 Y = ba_intervention, X = table_intervention_1_ba_2024)

#Now we do it for the breastmilk intervention samples

#restricting log feature table to only be the intervention baby samples
mb_intervention_id <- c("MB001B", "MB001E", "MB002B", "MB002E", "MB003B", "MB003E", "MB005B", "MB005E")
table_intervention_log_mb_2024 <- features_jamie_2024_log[mb_intervention_id, ]
View(table_intervention_log_mb_2024)
table_intervention_1_mb_2024 <- features_jamie_2024_1[mb_intervention_id, ]

#Clr function
clr_intervention_mb_2024 <- apply(table_intervention_log_mb_2024, 2, function(x) x- rowMeans(table_intervention_log_mb_2024))

#Determining the penalization parameter (lambda)
intervention_mb_2024.test_clrlasso <- glmnet(x = clr_intervention_mb_2024, y = ba_intervention, family = "binomial", nlambda = 50)
  # Warning message:
  #   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
  #               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_mb_2024.test_clrlasso, xvar = "lambda", label = T)
intervention_mb_2024.test_clrlasso
#we will choose a lambda of 0.007 (same as the lambda we used in the 2020 version)
#with 6 factors explaining 98% of the variation

intervention_mb_2024_clrlasso <- glmnet(x = clr_intervention_mb_2024, y = mom_bm_intervention, family = "binomial", lambda = 0.007)
  # Warning message:
  #   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
  #               one multinomial or binomial class has fewer than 8  observations; dangerous ground
intervention_mb_2024.results_clrlasso <- glmnet_wrapper(intervention_mb_2024_clrlasso, X = clr_intervention_mb_2024)

intervention_mb_2024.results_clrlasso$numVarSelect
  #6
intervention_mb_2024.results_clrlasso$varSelect
  # [1] "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__"                                
  # [2] "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae"                
  # [3] "k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__" 
  # [4] "k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__" 
  # [5] "k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Streptophyta; f__; g__; s__"                                      
  # [6] "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__"

#selbal like plot for clr-lasso
intervention_mb_2024.clr_pos <- intervention_mb_2024.results_clrlasso$posCoefSelect
intervention_mb_2024.clr_neg <- intervention_mb_2024.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_mb_2024.clr_pos), 
                 neg.names = names(intervention_mb_2024.clr_neg),
                 Y = mom_bm_intervention, X = table_intervention_1_mb_2024)


#Now we do it for the mom's feces intervention samples

#restricting log feature table to only be the intervention mom fecal samples
ma_intervention_id <- c("MA001B", "MA001E", "MA002B", "MA002E", "MA005B", "MA005E")
table_intervention_log_ma_2024 <- features_jamie_2024_log[ma_intervention_id, ]
View(table_intervention_log_ma_2024)
table_intervention_1_ma_2024 <- features_jamie_2024_1[ma_intervention_id, ]
View(table_intervention_1_ma_2024)

#Clr function
clr_intervention_ma_2024 <- apply(table_intervention_log_ma_2024, 2, function(x) x- rowMeans(table_intervention_log_ma_2024))
View(clr_intervention_ma_2024)

#Determining the penalization parameter (lambda)
intervention_ma_2024.test_clrlasso <- glmnet(x = clr_intervention_ma_2024, y = ma_intervention, family = "binomial", nlambda = 50)
# In plotCoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,  :
#               1 or less nonzero coefficients; glmnet plot is not meaningful
# Warning message:
#   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_ma_2024.test_clrlasso, xvar = "lambda", label = T)
intervention_ma_2024.test_clrlasso
#we will choose a lambda of 0.001
#with 2 factors explaining 99.9% of the variation

intervention_ma_2024_clrlasso <- glmnet(x = clr_intervention_ma_2024, y = ma_intervention, family = "binomial", lambda = 0.001)
# Warning message:
#   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#               one multinomial or binomial class has fewer than 8  observations; dangerous ground
intervention_ma_2024.results_clrlasso <- glmnet_wrapper(intervention_ma_2024_clrlasso, X = clr_intervention_ma_2024)

intervention_ma_2024.results_clrlasso$numVarSelect
#2
intervention_ma_2024.results_clrlasso$varSelect
# [1] "k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__"
# [2] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__"   

#selbal like plot for clr-lasso
intervention_ma_2024.clr_pos <- intervention_ma_2024.results_clrlasso$posCoefSelect
intervention_ma_2024.clr_neg <- intervention_ma_2024.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ma_2024.clr_pos), 
                 neg.names = names(intervention_ma_2024.clr_neg),
                 Y = ma_intervention, X = table_intervention_1_ma_2024)
