#Following the "Variable selection in microbiome compositional data analysis: tutorial"
#available at https://malucalle.github.io/Microbiome-Variable-Selection/index.html

system('git clone https://github.com/UVic-omics/CoDA-Penalized-Regression')

devtools::install_github(repo = "malucalle/selbal")
devtools::install_github(repo = "jbisanz/qiime2R")

library(knitr)
library(glmnet)
library(selbal)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(ggforce)
library(grid)
library(qiime2R)
library(tidyverse)

# source coda-lasso functions
source(file = './CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')
# build in functions
source(file = './CoDA-Penalized-Regression/R/functions.R')

#Importing metadata file, feature table, and taxonomy
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA")
metadata_mom_breastmilk <- read.table("meema_mom_bm_metadata.tsv", sep = "\t", stringsAsFactors = TRUE)
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA")
features_2020 <- read_qza("table-dada2.qza")
View(features_2020)
taxonomy_jamie_2020 <- read_qza("taxonomy_jamie_2020.qza")

#making taxonomy dictionary for CoDA outputs
taxonomy_jamie_dict <- setNames(taxonomy_jamie_2020$data$Taxon, taxonomy_jamie_2020$data$Feature.ID) 

#flipping rows/columns for feature data to match format in tutorial for selbal
inverted_table_2020 <- data.frame(t(features_2020$data))
View(inverted_table_2020)
write.csv(inverted_table_2020, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/inverted_table_2020.csv")

#making feature table with only the breastmilk samples from the 4 intervention individuals 
##for which we have both time points

intervention_bm_samples <- c("MB001B", "MB001E", "MB002B", "MB002E", "MB003B", "MB003E", "MB005B", "MB005E")
table_intervention_bm <- inverted_table_2020 %>% filter(row.names(inverted_table_2020) %in% intervention_bm_samples)
View(table_intervention_bm)

#making a y variable for the selbal analysis that gives the timepoints for the intervention samples
mom_bm_intervention <- factor(c("Baseline", "Endpoint", "Baseline", "Endpoint","Baseline", "Endpoint","Baseline", "Endpoint"))

##SELBAL##
selbal_bm_intervention <- selbal(x = table_intervention_bm, y = mom_bm_intervention, maxV = 12, 
                      logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)

taxonomy_jamie_dict["a6c0f31db841c236352967a9f4d5eed6"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__"
taxonomy_jamie_dict["6ec29d08a4d75314bb3c70989c1bb2d0"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__SMB53; s__"
taxonomy_jamie_dict["3a04dc65d8db8362d69f6a0fe36ceae5"]
  #"k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Nitrosomonadales; f__Nitrosomonadaceae; g__; s__"

#Running the rest of the selbal code
bm_intervention.results_selbal <- selbal_wrapper(result = selbal_bm_intervention, X = table_intervention_bm)
bm_intervention.selbal_pos <- bm_intervention.results_selbal$posVarSelect
bm_intervention.selbal_neg <- bm_intervention.results_selbal$negVarSelect
selbal_like_plot(pos.names = bm_intervention.selbal_pos, neg.names = bm_intervention.selbal_neg, 
                 Y = mom_bm_intervention, selbal = TRUE, 
                 FINAL.BAL = bm_intervention.results_selbal$finalBal)

##CLR-LASSO##

#Prepping the data by adding 1 and taking the log of all features
table_intervention_bm_1 <- table_intervention_bm + 1
table_intervention_bm_log <- log(table_intervention_bm_1)

#Clr function
clr_intervention_bm <- apply(table_intervention_bm_log, 2, function(x) x- rowMeans(table_intervention_bm_log))

#Determining the penalization parameter (lambda)
intervention_bm.test_clrlasso <- glmnet(x = clr_intervention_bm, y = mom_bm_intervention, family = "binomial", nlambda = 50)
# Warning message:
#   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_bm.test_clrlasso, xvar = "lambda", label = T)
intervention_bm.test_clrlasso
#a lambda of 0.004 explains 99% of the variation with 7 different factors; a lambda
#of 0.007 explains 98% of variation with 6 factors- but also a bigger lambda with
#the same number of factors explains less variation. 
#going to choose 0.06 as my lambda value
intervention_bm_clrlasso <- glmnet(x = clr_intervention_bm, y = mom_bm_intervention, family = "binomial", lambda = 0.007)
intervention_bm.results_clrlasso <- glmnet_wrapper(intervention_bm_clrlasso, X = clr_intervention_bm)
intervention_bm.results_clrlasso$numVarSelect
intervention_bm.results_clrlasso$varSelect
# [1] "X72c2f61f1aa3e438bfc3db3be87389ae" "X3da56361f77134e247c26b63e25fafa9" "a6c0f31db841c236352967a9f4d5eed6" 
# [4] "ee5f6e5409521f03c92ca643ffe3d187"  "X23b4285e223e34b174368d3cb87aced3" "e5363b0415557e975388806a6c12e8a4" 
taxonomy_jamie_dict["72c2f61f1aa3e438bfc3db3be87389ae"]
  #"k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__"
taxonomy_jamie_dict["3da56361f77134e247c26b63e25fafa9"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae"
taxonomy_jamie_dict["a6c0f31db841c236352967a9f4d5eed6"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__" 
taxonomy_jamie_dict["ee5f6e5409521f03c92ca643ffe3d187"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__" 
taxonomy_jamie_dict["23b4285e223e34b174368d3cb87aced3"]
  #"k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Streptophyta; f__; g__; s__" 
taxonomy_jamie_dict["e5363b0415557e975388806a6c12e8a4"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__" 

#selbal like plot for clr-lasso
intervention_bm.clr_pos <- intervention_bm.results_clrlasso$posCoefSelect
intervention_bm.clr_nef <- intervention_bm.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_bm.clr_pos), 
                 neg.names = names(intervention_bm.clr_nef),
                 Y = mom_bm_intervention, X = table_intervention_bm_1)


##CODA-LASSO##
#Use the table_1, not table_log because the function takes the log for you
lambdaRange_codalasso(X = table_intervention_bm_1, y = mom_bm_intervention, lambdaSeq = seq(0.7, 1, 0.01))
#lambda of 0.89 gives seven features and explains 100% of the variation
codalasso_intervention_bm_0.89 <- coda_logistic_lasso(X = table_intervention_bm_1, y = mom_bm_intervention, lambda = 0.89)
  #Warning messages:
  # 1: glm.fit: algorithm did not converge 
  # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
  # 3: glm.fit: algorithm did not converge 
  # 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
  # 5: glm.fit: fitted probabilities numerically 0 or 1 occurred 
  # 6: glm.fit: algorithm did not converge 
  # 7: glm.fit: fitted probabilities numerically 0 or 1 occurred 
intervention_bm.results_codalasso_0.89 <- coda_lasso_wrapper(result = codalasso_intervention_bm_0.89, X = table_intervention_bm_1)
intervention_bm.results_codalasso_0.89$numVarSelect
  #[1] 7
intervention_bm.results_codalasso_0.89$varSelect
  # [1] "a6c0f31db841c236352967a9f4d5eed6"  "aa977eb3ddebb896b1cb83ed13a501c5"  "X68ecf4cad684ba06c72dac6474479948"
  # [4] "X3da56361f77134e247c26b63e25fafa9" "efa65c4b52cfbdfc18fb223856de2c8b"  "d5ecb9dfd7960639bbd890dbec595043" 
  # [7] "X7d8be4a9e570c9db055c427977a98cbe"
taxonomy_jamie_dict["a6c0f31db841c236352967a9f4d5eed6"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Allobaculum; s__" 
taxonomy_jamie_dict["aa977eb3ddebb896b1cb83ed13a501c5"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__" 
taxonomy_jamie_dict["68ecf4cad684ba06c72dac6474479948"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__" 
taxonomy_jamie_dict["3da56361f77134e247c26b63e25fafa9"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae"
taxonomy_jamie_dict["efa65c4b52cfbdfc18fb223856de2c8b"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__"
taxonomy_jamie_dict["d5ecb9dfd7960639bbd890dbec595043"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__" 
taxonomy_jamie_dict["7d8be4a9e570c9db055c427977a98cbe"]
  #"k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae" 

#selbal like plot for coda-lasso
intervention_bm.coda_pos <- intervention_bm.results_codalasso_0.89$posCoefSelect
intervention_bm.coda_neg <- intervention_bm.results_codalasso_0.89$negCoefSelect
selbal_like_plot(pos.names = names(intervention_bm.coda_pos),
                 neg.names = names(intervention_bm.coda_neg),
                 Y = mom_bm_intervention, X = table_intervention_bm_1)


##Upset plot##
intervention_bm.select <- list(selbal = bm_intervention.results_selbal$varSelect,
                            clr_lasso = intervention_bm.results_clrlasso$varSelect,
                            coda_lasso = intervention_bm.results_codalasso_0.89$varSelect)
intervention_bm.select.upsetR <- fromList(intervention_bm.select)
upset(as.data.frame(intervention_bm.select.upsetR), main.bar.color = "gray36",
      sets.bar.color = color[c(1:2,5)], matrix.color = "gray36",
      order.by = "freq", empty.intersections = "on",
      queries = list(list(query = intersects, params = list("selbal"),
                          color = color[5], active = T),
                     list(query = intersects, params = list("clr_lasso"),
                          color = color[2], active = T),
                     list(query = intersects, params = list("coda_lasso"),
                          color = color[1], active = T)),
      text.scale = c(2, 2, 2, 2, 2, 2))

#only 1 feature is common among the 3 analyses is a6c0f31db841c236352967a9f4d5eed6, 
#which is allobaculum (genus)
