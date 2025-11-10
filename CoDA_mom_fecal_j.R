##I want to redo the CoDA analysis from the moms' fecal samples because I originally
##used Kiersten's taxonomy for the CoDA but redid it with mine when I did Sourcetracker
##because there was cross-version compatibility issues with the qiime2 files. 
##I'm worried that the important taxa identified between the different analyses is
##because her taxonomy is different than mine.

#I'm also only going to include the individuals for which we have both baseline and endpoint
#samples because that's what we used in the sourcetracker analysis

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

#uploading new taxonomy
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA")
taxonomy_jamie <- read_qza("taxonomy_jamie.qza")
taxonomy_jamie_2020 <- read_qza("taxonomy_jamie_2020.qza")

#making dictionary
taxonomy_jamie_dict <- setNames(taxonomy_jamie_2020$data$Taxon, taxonomy_jamie_2020$data$Feature.ID) 

#Importing feature table
features_2020 <- read_qza("table-dada2.qza")
View(features_2020$data)

#flipping rows/columns for feature data to match format in tutorial for selbal
inverted_table_2020 <- data.frame(t(features_2020$data))
View(inverted_table_2020)
write.csv(inverget ted_table_2020, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/inverted_table_2020.csv")
#The file looks good- column 1 is sample names and the rest are the features

#Now we're going to run selbal for the intervention group baseline vs. endpoint
#Making feature table
table_intervention_ma_2020 <- inverted_table_2020 %>% filter(row.names(inverted_table_2020) %in% 
                                                               c("MA001B", "MA001E", "MA002B", "MA002E", "MA003B", "MA003E", "MA005B", "MA005E"))
View(table_intervention_ma_2020)
#So we ran into an issue because sample MA003E wasn't included in the table, which
#didn't make sense. After looking back, there was no .fastq file for this sample so it
#hasn't been included in any analysis that I've run. 
##After looking into it more, I don't think this ever got sequenced, potentially because
##there wasn't enough of a fecal sample to begin with.

#Making y variable with timepoints as factor
row.names(table_intervention_ma_2020)
ma_intervention <- factor(c("Baseline", "Endpoint", "Baseline", "Endpoint",
                                   "Baseline", "Baseline", "Endpoint"))

selbal_intervention_ma <- selbal(x = table_intervention_ma_2020, y = ma_intervention, maxV = 12, 
                              logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)

#Running the rest of the selbal code
intervention_ma.results_selbal <- selbal_wrapper(result = selbal_intervention, X = table_intervention_ma_2020)
intervention_ma.selbal_pos <- intervention_ma.results_selbal$posVarSelect
intervention_ma.selbal_neg <- intervention_ma.results_selbal$negVarSelect
selbal_like_plot(pos.names = intervention_ma.selbal_pos, neg.names = intervention_ma.selbal_neg, 
                 Y = ma_intervention, selbal = TRUE, 
                 FINAL.BAL = intervention_ma.results_selbal$finalBal)

taxonomy_jamie_dict["7f8de42223057d36c757176022507598"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae; g__Gallibacterium; s__" 
taxonomy_jamie_dict["777066596250864f8f3c483a19583c8a"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales" 
taxonomy_jamie_dict["f2def463625868c5c9206b22ca0df43a"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__; s__" 

##CLR-LASSO##

#Prepping the data by adding 1 and taking the log of all features
table_intervention_ma_2020_1 <- table_intervention_ma_2020 + 1
table_intervention_ma_2020_log <- log(table_intervention_ma_2020_1)

#Clr function
clr_intervention_ma_2020 <- apply(table_intervention_ma_2020_log, 2, function(x) x- rowMeans(table_intervention_ma_2020_log))

#Determining the penalization parameter (lambda)
intervention_ma_2020.test_clrlasso <- glmnet(x = clr_intervention_ma_2020, y = ma_intervention, family = "binomial", nlambda = 100)
# Warning message:
#   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_ma_2020.test_clrlasso, xvar = "lambda", label = T)
intervention_ma_2020.test_clrlasso
#a lambda of 0.00651 explains 99.0% of the variation with 1  factor- have to choose a much smaller
#lambda to get more than 1 factor
intervention_ma_2020_clrlasso <- glmnet(x = clr_intervention_ma_2020, y = ma_intervention, family = "binomial", lambda = 0.0001)
#Warking: one multinomial or binomial class has fewer than 8 observations; dangerous ground
intervention_ma_2020.results_clrlasso <- glmnet_wrapper(intervention_ma_2020_clrlasso, X = clr_intervention_ma_2020)
intervention_ma_2020.results_clrlasso$numVarSelect
intervention_ma_2020.results_clrlasso$varSelect
  #[1] "f2def463625868c5c9206b22ca0df43a"  "X6195e5ef381b20047179f7d4dcff7952" "X0cb6ba500c97a4dd843b414c972739d6"
taxonomy_jamie_dict["f2def463625868c5c9206b22ca0df43a"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__; s__" 
taxonomy_jamie_dict["6195e5ef381b20047179f7d4dcff7952"]
  #"k__Bacteria"
taxonomy_jamie_dict["0cb6ba500c97a4dd843b414c972739d6"]
  #"k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__longum"

#selbal like plot for clr-lasso
intervention_ma_2020.clr_pos <- intervention_ma_2020.results_clrlasso$posCoefSelect
intervention_ma_2020.clr_nef <- intervention_ma_2020.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ma_2020.clr_pos), 
                 neg.names = names(intervention_ma_2020.clr_nef),
                 Y = ma_intervention, X = table_intervention_ma_2020_1)


##CODA-LASSO##
#Use the table_1, not table_log because the function takes the log for you
lambdaRange_codalasso(X = table_intervention_ma_2020_1, y = ma_intervention, lambdaSeq = seq(1, 2, 0.05))
#lambda of 1.35 gives five features and explains 100% of the variation
codalasso_intervention_ma_2020 <- coda_logistic_lasso(X = table_intervention_ma_2020_1, y = ma_intervention, lambda = 1.35)
#Warning messages:
# 1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 3: glm.fit: algorithm did not converge 
# 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 5: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 6: glm.fit: algorithm did not converge 
# 7: glm.fit: fitted probabilities numerically 0 or 1 occurred 
intervention_ma_2020.results_codalasso <- coda_lasso_wrapper(result = codalasso_intervention_ma_2020, X = table_intervention_ma_2020_1)
intervention_ma_2020.results_codalasso$numVarSelect
#[1] 5
intervention_ma_2020.results_codalasso$varSelect
#[1] "f2def463625868c5c9206b22ca0df43a"  "X691fcc90248f12d2304e3376c38fba11" "X849c75ff46a2f0c7ae74eac904fa9752"
#[4] "X6fceba9d0045ef080635acee41fd5965" "X9d79ae842e336b7646b18760156f6e37"
taxonomy_jamie_dict["f2def463625868c5c9206b22ca0df43a"]
  #"k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__; s__" 
taxonomy_jamie_dict["691fcc90248f12d2304e3376c38fba11"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__[Tissierellaceae]; g__Anaerococcus; s__"
taxonomy_jamie_dict["849c75ff46a2f0c7ae74eac904fa9752"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae" 
taxonomy_jamie_dict["6fceba9d0045ef080635acee41fd5965"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__"
taxonomy_jamie_dict["9d79ae842e336b7646b18760156f6e37"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Epulopiscium; s__" 


#selbal like plot for coda-lasso
intervention_ma_2020.coda_pos <- intervention_ma_2020.results_codalasso$posCoefSelect
intervention_ma_2020.coda_neg <- intervention_ma_2020.results_codalasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ma_2020.coda_pos),
                 neg.names = names(intervention_ma_2020.coda_neg),
                 Y = ma_intervention, X = table_intervention_ma_2020_1)


##Upset plot##
intervention_ma_2020.select <- list(selbal = intervention_ma.results_selbal$varSelect,
                               clr_lasso = intervention_ma_2020.results_clrlasso$varSelect,
                               coda_lasso = intervention_ma_2020.results_codalasso$varSelect)
intervention_ma_2020.select.upsetR <- fromList(intervention_ma_2020.select)
upset(as.data.frame(intervention_ma_2020.select.upsetR), main.bar.color = "gray36",
      sets.bar.color = color[c(1:2,5)], matrix.color = "gray36",
      order.by = "freq", empty.intersections = "on",
      queries = list(list(query = intersects, params = list("selbal"),
                          color = color[5], active = T),
                     list(query = intersects, params = list("clr_lasso"),
                          color = color[2], active = T),
                     list(query = intersects, params = list("coda_lasso"),
                          color = color[1], active = T)),
      text.scale = c(2, 2, 2, 2, 2, 2))




