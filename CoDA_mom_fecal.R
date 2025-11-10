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

# source coda-lasso functions
source(file = './CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')
# build in functions
source(file = './CoDA-Penalized-Regression/R/functions.R')

#Importing metadata files
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA-Fastq")
metadata <- read.csv("meema_metadata.csv", stringsAsFactors = TRUE)
metadata_moms <- read.csv("meema_metadata_mom.csv", stringsAsFactors = TRUE)
setwd("//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/Kiersten")
metadata_mom_fecal <- read.table("meema_mom_fecal_metadata.tsv", sep = "\t", stringsAsFactors = TRUE)
metadata_mom_breastmilk <- read.table("meema_mom_bm_metadata.tsv", sep = "\t", stringsAsFactors = TRUE)
metadata_baby_fecal <- read.table("meema_baby_fecal_metadata.tsv", sep = "\t", stringsAsFactors = TRUE)

#Importing feature table and taxonomy for the fecal samples
feature_table_mom_fec <- read_qza("MEEMA_M_Fec_table.qza")
taxonomy_mom_fec <- read_qza("MEEMA_M_Fec_taxonomy.qza")

#flipping rows/columns for feature data to match format in tutorial for selbal
inverted_table <- data.frame(t(feature_table_mom_fec$data))
colnames(inverted_table) <- feature_table_mom_fec$data[, 1]
head(inverted_table)
summary(inverted_table)
write.csv(inverted_table, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/inverted_table.csv")
#The file looks good- column 1 is sample names and the rest are the features
#This is only from the the moms' fecal samples though don't forget

#making feature table for only the control samples at the different timepoints
rownames(inverted_table)
library(dplyr)
table_ctrl <- inverted_table %>% filter(row.names(inverted_table) %in% c("MA004B", "MA004E", "MA007B"))
row.names(table_ctrl)

#making dataframe with the sample names and the timepoints
mom_fecal_ctrl <- metadata_mom_fecal %>% filter(metadata_mom_fecal$V1 %in% c("MA004B", "MA004E", "MA007B"))
mom_fecal_ctrl <- subset(mom_fecal_ctrl, select = -c(V2, V3, V4))
mom_fecal_ctrl

#now actually running selbel for the control at baseline and endpoint
selbal_ctrl <- selbal(x = table_ctrl, y = mom_fecal_ctrl, maxV = 12, 
                       logit.acc = 'Dev', draw = T, logt = TRUE)
#This threw an error: "GBM method: not enough information to compute t hyper-parameter"

#So trying to remove columns with 0 for all samples
table_ctrl <- Filter(function(x) !all(x == 0), table_ctrl)

selbal_ctrl <- selbal(x = table_ctrl, y = mom_fecal_ctrl, maxV = 12, 
                      logit.acc = 'Dev', draw = T, logt = TRUE)
#Same error thrown, but removed the warnings so I guess that's a good thing.
# Error in cmultRepl(x, z.delete = FALSE, suppress.print = T) : 
#   GBM method: not enough information to compute t hyper-parameter,
# probably there are columns with < 2 positive values.

#I'm going to mess with the zero.rep function to see if that fixes the error in BM calcs
selbal_ctrl <- selbal(x = table_ctrl, y = mom_fecal_ctrl, maxV = 12, 
                      logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)
# Error in data.frame(..., check.names = FALSE) : 
#   arguments imply differing number of rows: 3, 0
summary(mom_fecal_ctrl)
#I think I need a factor for y instead of a df
row.names(mom_fecal_ctrl) <- c("MA004B", "MA004E", "MA007B")
row.names(mom_fecal_ctrl)
mom_fecal_ctrl <- factor(c("Baseline", "Endpoint", "Baseline"))

selbal_ctrl <- selbal(x = table_ctrl, y = mom_fecal_ctrl, maxV = 12, 
                      logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)
                
#This is good, but I don't want the feature names there- I want to have the genus or species for the feature.

#Removes all taxonomy classes before final semicolon, thus keeping the most specific known info.
short_tax_mom_fec <- sub(".*;(.*)", "\\1", taxonomy_mom_fec$data$Taxon)

#Making dictionary with feature and taxonomy
mom_fec_dict <- setNames(taxonomy_mom_fec$data$Taxon, taxonomy_mom_fec$data$Feature.ID) 

#Making dictionary with feature and short taxonomy
mom_fec_dict_short <- setNames(short_tax_mom_fec, taxonomy_mom_fec$data$Feature.ID)

#Running the rest of the selbal code
ctrl.results_selbal <- selbal_wrapper(result = selbal_ctrl, X = table_ctrl)
ctrl.selbal_pos <- ctrl.results_selbal$posVarSelect
ctrl.selbal_neg <- ctrl.results_selbal$negVarSelect
selbal_like_plot(pos.names = ctrl.selbal_pos, neg.names = ctrl.selbal_neg, 
                 Y = mom_fecal_ctrl, selbal = TRUE, 
                 FINAL.BAL = ctrl.results_selbal$finalBal)
# Warning messages:
#   1: Groups with fewer than two data points have been dropped. 
# 2: In max(ids, na.rm = TRUE) :
#   no non-missing arguments to max; returning -Inf

mom_fec_dict_short[ctrl.results_selbal$posVarSelect]
#" g__Lactobacillus" 
mom_fec_dict_short[ctrl.results_selbal$negVarSelect]
#" s__CAG-313 sp000433035"

#Now we're going to run selbal for the intervention group baseline vs. endpoint
#Making feature table
table_intervention <- inverted_table[!(row.names(inverted_table) %in% c("MA004B", "MA004E", "MA007B")), ]
#Making y variable with timepoints as factor
row.names(table_intervention)
mom_fecal_intervention <- factor(c("Baseline", "Endpoint", "Baseline", "Endpoint",
                                  "Baseline", "Baseline", "Endpoint", "Baseline", "Baseline"))
#Removing  columns with 0 for all samples
table_intervention <- Filter(function(x) !all(x == 0), table_intervention)

selbal_intervention <- selbal(x = table_intervention, y = mom_fecal_intervention, maxV = 12, 
                      logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)
#Running the rest of the selbal code
intervention.results_selbal <- selbal_wrapper(result = selbal_intervention, X = table_intervention)
intervention.selbal_pos <- intervention.results_selbal$posVarSelect
intervention.selbal_neg <- intervention.results_selbal$negVarSelect
selbal_like_plot(pos.names = intervention.selbal_pos, neg.names = intervention.selbal_neg, 
                 Y = mom_fecal_intervention, selbal = TRUE, 
                 FINAL.BAL = intervention.results_selbal$finalBal)

mom_fec_dict["ddf8c967a882b6a5d50cf4039f7b7966"]
# "d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Oscillospirales; f__Acutalibacteraceae; g__UMGS1071; s__UMGS1071 sp900548305"

mom_fec_dict["3979b271c664505fb0cbd44c6900617a"]
# "d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Copromonas; s__Copromonas sp900066535" 
mom_fec_dict["627a28ae2ce79a8efbdfbd712528b1f2"]
# "d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Oscillospirales; f__Ruminococcaceae; g__Anaerofilum_74150; s__"

#Now selbal for the control vs. intervention at baseline
#Making feature table
table_baseline <- inverted_table %>% filter(row.names(inverted_table) %in% c("MA001B", "MA002B", "MA003B", "MA004B", "MA005B", "MA006B", "MA007B", "MA008B"))
#Making y variable with control vs. intervention as factor
mom_fecal_baseline <- factor(c("I", "I", "I", "C", "I", "I", "C", "I"))
#Removing  columns with 0 for all samples
table_baseline <- Filter(function(x) !all(x == 0), table_baseline)                                                                            

#running selbal function
selbal_baseline <- selbal(x = table_baseline, y = mom_fecal_baseline, maxV = 12, 
                              logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)

#Running the rest of the selbal code
baseline.results_selbal <- selbal_wrapper(result = selbal_baseline, X = table_baseline)
baseline.selbal_pos <- baseline.results_selbal$posVarSelect
baseline.selbal_neg <- baseline.results_selbal$negVarSelect
selbal_like_plot(pos.names = baseline.selbal_pos, neg.names = baseline.selbal_neg, 
                 Y = mom_fecal_baseline, selbal = TRUE, 
                 FINAL.BAL = baseline.results_selbal$finalBal)

#Neg var
mom_fec_dict["8da2eff472c19c4828f0d73dc3933c5e"]
# "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides_H; s__Bacteroides_H cellulosilyticus" 

#Pos var
mom_fec_dict["b37d168846a2f72adabfc8f0edcb7d55"]
# "d__Bacteria; p__Firmicutes_D; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Limosilactobacillus; s__Limosilactobacillus fermentum" 

#Now selbal for the control vs. intervention at endpoint
#Making feature table
table_endpoint <- inverted_table %>% filter(row.names(inverted_table) %in% c("MA001E", "MA002E", "MA004E", "MA005E"))
#Making y variable with control vs. intervention as factor
mom_fecal_endpoint <- factor(c("I", "I", "C", "I"))
#Removing  columns with 0 for all samples
table_endpoint <- Filter(function(x) !all(x == 0), table_endpoint)                                                                            

#running selbal function
selbal_endpoint <- selbal(x = table_endpoint, y = mom_fecal_endpoint, maxV = 12, 
                          logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)

#Running the rest of the selbal code
endpoint.results_selbal <- selbal_wrapper(result = selbal_endpoint, X = table_endpoint)
endpoint.selbal_pos <- endpoint.results_selbal$posVarSelect
endpoint.selbal_neg <- endpoint.results_selbal$negVarSelect
selbal_like_plot(pos.names = endpoint.selbal_pos, neg.names = endpoint.selbal_neg, 
                 Y = mom_fecal_endpoint, selbal = TRUE, 
                 FINAL.BAL = endpoint.results_selbal$finalBal)

#Neg var
mom_fec_dict["37f9d87b0455de5aee8807f90aec4e92"]
# "d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae" 

#Pos var
mom_fec_dict["bd2dbb0c01a0db803b60bd021f1a540b"]
# "d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Eubacteriales_258471; f__Eubacteriaceae_258449; g__Eubacterium_O_258270" 


# The next part requires non-zero values for counts, so I'm thinking maybe I should try that for selbal and see if I get 
# any different outcomes-- it didn't
table_intervention_1 <- table_intervention + 1

selbal_intervention_1 <- selbal(x = table_intervention_1, y = mom_fecal_intervention, maxV = 12, 
                              logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)
#Running the rest of the selbal code
intervention_1.results_selbal <- selbal_wrapper(result = selbal_intervention_1, X = table_intervention_1)
intervention_1.selbal_pos <- intervention_1.results_selbal?$posVarSelect
intervention_1.selbal_neg <- intervention_1.results_selbal$negVarSelect
selbal_like_plot(pos.names = intervention_1.selbal_pos, neg.names = intervention_1.selbal_neg, 
                 Y = mom_fecal_intervention, selbal = TRUE, 
                 FINAL.BAL = intervention_1.results_selbal$finalBal)
##This didn't change the plot, so wouldn't help to redo for all the analysis 


###Clr-lasso time###
#Requires at least 1 for all quantities
table_baseline_1 <- table_baseline + 1
table_ctrl_1 <- table_ctrl + 1
table_endpoint_1 <- table_endpoint + 1
table_intervention_1 <- table_intervention + 1

#now taking the log of all tables
table_baseline_log <- log(table_baseline_1)
table_endpoint_log <- log(table_endpoint_1)
table_ctrl_log <- log(table_ctrl_1)
table_intervention_log <- log(table_intervention_1)

#clr function for each of the log tables
clr_baseline <- apply(table_baseline_log, 2, function(x) x - rowMeans(table_baseline_log))
clr_endpoint <- apply(table_endpoint_log, 2, function(x) x - rowMeans(table_endpoint_log))
clr_ctrl <- apply(table_ctrl_log, 2, function(x) x - rowMeans(table_ctrl_log))
clr_intervention <- apply(table_intervention_log, 2, function(x) x- rowMeans(table_intervention_log))

# y needs to be numeric for this analysis, so I have to change the y tables into numeric values
##Baseline = 0 ; Endpoint = 1 ; Control = 10 ; Intervention = 11
mom_fecal_baseline_num <- c(11, 11, 11, 10, 11, 11, 10, 11)
mom_fecal_endpoint_num <- c(11, 11, 10, 11)
mom_fecal_ctrl_num <- c(0, 1, 0)
mom_fecal_intervention_num <- c(0, 1, 0, 1, 0, 0, 1, 0, 0)

#Determining the penalization parameter (lambda)
baseline.test_clrlasso <- glmnet(x = clr_baseline, y = mom_fecal_baseline_num, family = "binomial", nlambda = 30)
#"Warning message:
#In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#  one multinomial or binomial class has fewer than 8  observations; dangerous ground"
plot(baseline.test_clrlasso, xvar = "lambda", label = T)
baseline.test_clrlasso
#changing the nlambda from 30 -> 100 didn't change the overall percent variation explained, so 
#I will keep it at 30
#The tutorial says that they set the optimal discrimination to be 12 variables because that's how
#many they found in the selbal. We only saw 2-4, so I'm going to choose between those depending on the dataset.
#The optimal lambda value for the baseline test is about 0.00972, which I'm just going to round to 0.009.
baseline_clrlasso <- glmnet(x = clr_baseline, y = mom_fecal_baseline_num, family = "binomial", lambda = 0.009)
baseline.results_clrlasso <- glmnet_wrapper(baseline_clrlasso, X = clr_baseline)
baseline.results_clrlasso$numVarSelect
baseline.results_clrlasso$varSelect
# [1] "d8a37c9ad1d0f5e42c521eb6caeb03b4"  "X499a788d909315d1ae4628ec26d260c7"
# [3] "X82305c50e64befde41fa0ad19bb2977c" "X122ec7888a99827e071355d08c7fc4e5"
mom_fec_dict_short["d8a37c9ad1d0f5e42c521eb6caeb03b4"]
#g__Lactobacillus
mom_fec_dict_short["499a788d909315d1ae4628ec26d260c7"]
#s__Roseburia inulinivorans
mom_fec_dict_short["82305c50e64befde41fa0ad19bb2977c"]
#s__Akkermansia sp001580195
mom_fec_dict_short["122ec7888a99827e071355d08c7fc4e5"]
#" g__Staphylococcus"

#Everything looks good, so we are going to move on and run the workflow with the rest of the data
#Determining the penalization parameter (lambda)
endpoint.test_clrlasso <- glmnet(x = clr_endpoint, y = mom_fecal_endpoint_num, family = "binomial", nlambda = 30)
#can't do this with the endpoint because we only have one control sample from the endpoint.
#I'm not even going to try with the control because I know that we will run into the same problem
intervention.test_clrlasso <- glmnet(x = clr_intervention, y = mom_fecal_intervention_num, family = "binomial", nlambda = 30)
plot(intervention.test_clrlasso, xvar = "lambda", label = T)
intervention.test_clrlasso
#We get 6 features, which correlates to 99% deviation (idk how to interpret this yet)
#lambda = 0.00465
intervention_clrlasso <- glmnet(x = clr_intervention, y = mom_fecal_intervention_num,
                                  family = "binomial", lambda = 0.00465)
intervention.results_clrlasso <- glmnet_wrapper(intervention_clrlasso, X = clr_intervention)
intervention.results_clrlasso$numVarSelect
# [1] 6
intervention.results_clrlasso$varSelect
# [1] "X8000f3da896603a9f7346fa0f3dea39b" "X99ada31d7e865fd84371f5919f63f2f3"
# [3] "ddf8c967a882b6a5d50cf4039f7b7966"  "X17d5a20d7a14b42b9e92b10be66e6a48"
# [5] "X0ce9f8b3537b977bb3b63659d5bf08ad" "X38c1288a8a179b5dfd9cffb9d351232b"

mom_fec_dict_short["8000f3da896603a9f7346fa0f3dea39b"]
#s__AF33-28 sp003477885
mom_fec_dict_short["99ada31d7e865fd84371f5919f63f2f3"]
#s__Phocea massiliensis
mom_fec_dict_short["ddf8c967a882b6a5d50cf4039f7b7966"]
#s__UMGS1071 sp900548305
mom_fec_dict_short["17d5a20d7a14b42b9e92b10be66e6a48"]
#s__Lachnospira eligens
mom_fec_dict_short["0ce9f8b3537b977bb3b63659d5bf08ad"]
#c__Bacilli
mom_fec_dict_short["38c1288a8a179b5dfd9cffb9d351232b"]
#g__Veillonella_A


###Moving onto coda-lasso###
#Use the table_1, not table_log because the function takes the log for you
##baseline
lambdaRange_codalasso(X = table_baseline_1, y = mom_fecal_baseline, lambdaSeq = seq(1, 2, 0.1))
#Lambda = 1.3 gives us 4 variables selected; 1.0 gives 13 but the probability doesn't change from 2 to infitity
#So picking lambda for 2 num.selected, which is 1.5
codalasso_baseline <- coda_logistic_lasso(X = table_baseline_1, y = mom_fecal_baseline, lambda = 1.5)
# Warning messages:
#   1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 3: glm.fit: algorithm did not converge 
# 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 5: glm.fit: fitted probabilities numerically 0 or 1 occurred 

baseline.results_codalasso <- coda_lasso_wrapper(result = codalasso_baseline, X = table_baseline_1)
baseline.results_codalasso$numVarSelect
#2
baseline.results_codalasso$varSelect
# [1] "X82305c50e64befde41fa0ad19bb2977c" "X41a81e46174b2eef0a99567fc631c798"

mom_fec_dict_short["82305c50e64befde41fa0ad19bb2977c"]
#s__Akkermansia sp001580195
mom_fec_dict_short["41a81e46174b2eef0a99567fc631c798"]
#s__Phascolarctobacterium_A faecium

##endpoint
lambdaRange_codalasso(X = table_endpoint_1, y = mom_fecal_endpoint, lambdaSeq = seq(1, 2, 0.1))
#1.9 gives 3 and explains 100% of the deviation
codalasso_endpoint <- coda_logistic_lasso(X = table_endpoint_1, y = mom_fecal_endpoint, lambda = 1.9)
#Surprisingly no warnings
endpoint.results_codalasso <- coda_lasso_wrapper(result = codalasso_endpoint, X = table_endpoint_1)
endpoint.results_codalasso$numVarSelect
#3
endpoint.results_codalasso$varSelect
# [1] "ff62d2ad7a289fe535520e2aa8caf237"  "X25d27528ff1a7f36f8e26878a76c9c62"
# [3] "aa1b2718ac04d75614f022e5b777fa85" 

mom_fec_dict_short["ff62d2ad7a289fe535520e2aa8caf237"]
#s__Mediterraneibacter_A_155507 faecis
mom_fec_dict_short["25d27528ff1a7f36f8e26878a76c9c62"]
#s__CAG-83 sp000431575
mom_fec_dict_short["aa1b2718ac04d75614f022e5b777fa85"]
#g__Turicibacter

##control
lambdaRange_codalasso(X = table_ctrl_1, y = mom_fecal_ctrl, lambdaSeq = seq(1, 4, 0.1))
#lambda = 3.1 gives 3 features and explains 100% of the variation
codalasso_ctrl <- coda_logistic_lasso(X = table_ctrl_1, y = mom_fecal_ctrl, lambda = 3.1)
#No warnings for this one either
ctrl.results_codalasso <- coda_lasso_wrapper(result = codalasso_ctrl, X = table_ctrl_1)
ctrl.results_codalasso$numVarSelect
#3
ctrl.results_codalasso$varSelect
# [1] "X82305c50e64befde41fa0ad19bb2977c" "X48ddab680f67c18dcb89e9fb209896d4"
# [3] "c84815279a0c037d84a5dcec13c5020c" 

mom_fec_dict_short["82305c50e64befde41fa0ad19bb2977c"]
#s__Akkermansia sp001580195
mom_fec_dict_short["48ddab680f67c18dcb89e9fb209896d4"]
#s__Akkermansia muciniphila_D_776786
mom_fec_dict_short["c84815279a0c037d84a5dcec13c5020c"]
#g__Clostridium_T

##intervention
lambdaRange_codalasso(X = table_intervention_1, y = mom_fecal_intervention, lambdaSeq = seq(1, 2, 0.1))
#Lambda of 1.1 gives 6 features and explains 100% of the deviation
codalasso_intervention_1.1 <- coda_logistic_lasso(X = table_intervention_1, y = mom_fecal_intervention, lambda = 1.1)
# Warning messages:
#   1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 3: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
codalasso_intervention_1.2 <- coda_logistic_lasso(X = table_intervention_1, y = mom_fecal_intervention, lambda = 1.2)
#When I run it with lambda = 1.2 (4 features and only 50% deviation explained), I run into no errors
intervention.results_codalasso_1.1 <- coda_lasso_wrapper(result = codalasso_intervention_1.1, X = table_intervention_1)
intervention.results_codalasso_1.1$numVarSelect
#6
intervention.results_codalasso_1.1$varSelect
# [1] "X8b3d2e6fbcf39e8baf2f769af8698963" "X31c0dcde74a37e939116f4ee14b63d3a"
# [3] "de2610d3a69d945db4c0e8197cc47cad"  "X37f9d87b0455de5aee8807f90aec4e92"
# [5] "ffbf0d3f91f68ddcb32f6e20454d698e"  "X3979b271c664505fb0cbd44c6900617a"

mom_fec_dict_short["8b3d2e6fbcf39e8baf2f769af8698963"]
#s__Dialister invisus
mom_fec_dict_short["31c0dcde74a37e939116f4ee14b63d3a"]
#s__Ruminococcus_E sp003438075
mom_fec_dict_short["de2610d3a69d945db4c0e8197cc47cad"]
#s__Clostridium_Q_135853 saccharolyticum_A
mom_fec_dict_short["37f9d87b0455de5aee8807f90aec4e92"]
#f__Lachnospiraceae
mom_fec_dict_short["ffbf0d3f91f68ddcb32f6e20454d698e"]
#s__Lachnospira eligens
mom_fec_dict_short["3979b271c664505fb0cbd44c6900617a"]
#s__Copromonas sp900066535

intervention.results_codalasso_1.2 <- coda_lasso_wrapper(result = codalasso_intervention_1.2, X = table_intervention_1)
intervention.results_codalasso_1.2$numVarSelect
#4
intervention.results_codalasso_1.2$varSelect
# [1] "X31c0dcde74a37e939116f4ee14b63d3a" "X8b3d2e6fbcf39e8baf2f769af8698963"
# [3] "de2610d3a69d945db4c0e8197cc47cad"  "ffbf0d3f91f68ddcb32f6e20454d698e"
#4/6 from above: s__Dialister invisus, s__Ruminococcus_E sp003438075, 
#s__Clostridium_Q_135853 saccharolyticum_A, and s__Lachnospira eligens


##Going to make some plots to compare the results for all of the different methods
##for the intervention group because that is waht we are most interested in
#Upset plot
intervention.select <- list(selbal = intervention.results_selbal$varSelect,
                            clr_lasso = intervention.results_clrlasso$varSelect,
                            coda_lasso = intervention.results_codalasso_1.1$varSelect)
intervention.select.upsetR <- fromList(intervention.select)
upset(as.data.frame(intervention.select.upsetR), main.bar.color = "gray36",
      sets.bar.color = color[c(1:2,5)], matrix.color = "gray36",
      order.by = "freq", empty.intersections = "on",
      queries = list(list(query = intersects, params = list("selbal"),
                          color = color[5], active = T),
                     list(query = intersects, params = list("clr_lasso"),
                          color = color[2], active = T),
                     list(query = intersects, params = list("coda_lasso"),
                          color = color[1], active = T)),
     text.scale = c(2, 2, 2, 2, 2, 2))

#selbal like plot for clr-lasso
intervention.clr_pos <- intervention.results_clrlasso$posCoefSelect
intervention.clr_nef <- intervention.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention.clr_pos), 
                 neg.names = names(intervention.clr_nef),
                 Y = mom_fecal_intervention, X = table_intervention_1)

#selbal like plot for coda-lasso
intervention.coda_pos <- intervention.results_codalasso_1.1$posCoefSelect
intervention.coda_neg <- intervention.results_codalasso_1.1$negCoefSelect
selbal_like_plot(pos.names = names(intervention.coda_pos),
                 neg.names = names(intervention.coda_neg),
                 Y = mom_fecal_intervention, X = table_intervention_1)
