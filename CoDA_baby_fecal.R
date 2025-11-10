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
write.csv(inverted_table_2020, "//wsl.localhost/Ubuntu/home/jamiekaiser13/MEEMA/inverted_table_2020.csv")

#Now we're going to run selbal for the intervention group baseline vs. endpoint
#Making feature table with only baby fecal samples
table_intervention_ba_2020 <- inverted_table_2020 %>% filter(row.names(inverted_table_2020) %in% 
                                                               c("BA001B", "BA001E", "BA002B", "BA002E", "BA003B", "BA003E", "BA005B", "BA005E"))
View(table_intervention_ba_2020)
  #looks good

#Making y variable with timepoints as factor
row.names(table_intervention_ba_2020)
ba_intervention <- factor(c("Baseline", "Endpoint", "Baseline", "Endpoint",
                            "Baseline", "Endpoint", "Baseline", "Endpoint"))

selbal_intervention_ba <- selbal(x = table_intervention_ba_2020, y = ba_intervention, maxV = 12, 
                                 logit.acc = 'Dev', zero.rep = "one", draw = T, logt = TRUE)

#Running the rest of the selbal code
intervention_ba.results_selbal <- selbal_wrapper(result = selbal_intervention_ba, X = table_intervention_ba_2020)
intervention_ba.selbal_pos <- intervention_ba.results_selbal$posVarSelect
intervention_ba.selbal_neg <- intervention_ba.results_selbal$negVarSelect
selbal_like_plot(pos.names = intervention_ba.selbal_pos, neg.names = intervention_ba.selbal_neg, 
                 Y = ba_intervention, selbal = TRUE, 
                 FINAL.BAL = intervention_ba.results_selbal$finalBal)

taxonomy_jamie_dict["d1af5a4f25d0d352c615778dd8b42733"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__" 
taxonomy_jamie_dict["fa6297cc690bda4d77ddad4b720bd830"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Dialister; s__"
taxonomy_jamie_dict["ec30c0889697dc9b180d9c3aa13b430b"]
  #"k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__uniformis" 
taxonomy_jamie_dict["3c073cb98e82f33bfbe1e586326f062c"]
  #"k__Bacteria; p__Actinobacteria; c__Coriobacteriia; o__Coriobacteriales; f__Coriobacteriaceae; g__; s__"

##CLR-LASSO##

#Prepping the data by adding 1 and taking the log of all features
table_intervention_ba_2020_1 <- table_intervention_ba_2020 + 1
table_intervention_ba_2020_log <- log(table_intervention_ba_2020_1)

#Clr function
clr_intervention_ba_2020 <- apply(table_intervention_ba_2020_log, 2, function(x) x- rowMeans(table_intervention_ba_2020_log))

#Determining the penalization parameter (lambda)
intervention_ba_2020.test_clrlasso <- glmnet(x = clr_intervention_ba_2020, y = ba_intervention, family = "binomial", nlambda = 100)
# Warning message:
#   In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#               one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(intervention_ba_2020.test_clrlasso, xvar = "lambda", label = T)
intervention_ba_2020.test_clrlasso
#a lambda of 0.03152 explains 92.48% of the variation with 5 factors
intervention_ba_2020_clrlasso <- glmnet(x = clr_intervention_ba_2020, y = ba_intervention, family = "binomial", lambda = 0.03152)
#Warning: one multinomial or binomial class has fewer than 8 observations; dangerous ground
intervention_ba_2020.results_clrlasso <- glmnet_wrapper(intervention_ba_2020_clrlasso, X = clr_intervention_ba_2020)
intervention_ba_2020.results_clrlasso$numVarSelect
intervention_ba_2020.results_clrlasso$varSelect
  # [1] "d1af5a4f25d0d352c615778dd8b42733"  "X370c9b25114a189093b30487f62855a5" "a1842f623f2b17afd3460def32ead097" 
  # [4] "X28d93e6bfa906ffbd0da011b112ae983" "a4b9fa5b631a20f74556c49ee96bfe67" 
taxonomy_jamie_dict["d1af5a4f25d0d352c615778dd8b42733"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__" 
taxonomy_jamie_dict["370c9b25114a189093b30487f62855a5"]
  #"k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__europaeus" 
taxonomy_jamie_dict["a1842f623f2b17afd3460def32ead097"]
  #"k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__uniformis" 
taxonomy_jamie_dict["28d93e6bfa906ffbd0da011b112ae983"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__[Exiguobacteraceae]; g__; s__" 
taxonomy_jamie_dict["a4b9fa5b631a20f74556c49ee96bfe67"]
  #"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Dialister; s__" 

#selbal like plot for clr-lasso
intervention_ba_2020.clr_pos <- intervention_ba_2020.results_clrlasso$posCoefSelect
intervention_ba_2020.clr_nef <- intervention_ba_2020.results_clrlasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ba_2020.clr_pos), 
                 neg.names = names(intervention_ba_2020.clr_nef),
                 Y = ba_intervention, X = table_intervention_ba_2020_1)


##CODA-LASSO##
#Use the table_1, not table_log because the function takes the log for you
lambdaRange_codalasso(X = table_intervention_ba_2020_1, y = ba_intervention, lambdaSeq = seq(1, 2, 0.05))
#lambda of 1.55 gives four features and explains 100% of the variation
codalasso_intervention_ba_2020 <- coda_logistic_lasso(X = table_intervention_ba_2020_1, y = ba_intervention, lambda = 1.55)
#Warning messages:
# 1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 3: glm.fit: algorithm did not converge 
# 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 5: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 6: glm.fit: algorithm did not converge 
# 7: glm.fit: fitted probabilities numerically 0 or 1 occurred 
intervention_ba_2020.results_codalasso <- coda_lasso_wrapper(result = codalasso_intervention_ba_2020, X = table_intervention_ba_2020_1)
intervention_ba_2020.results_codalasso$numVarSelect
#[1] 4
intervention_ba_2020.results_codalasso$varSelect
  # [1] "d1af5a4f25d0d352c615778dd8b42733"  "X28d93e6bfa906ffbd0da011b112ae983" "b81ca62dc1c0e44c6f515c37a227d063" 
  # [4] "c363a0062f4a061d01b9c52da4011aa1" 
taxonomy_jamie_dict["d1af5a4f25d0d352c615778dd8b42733"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__" 
taxonomy_jamie_dict["28d93e6bfa906ffbd0da011b112ae983"]
  #"k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__[Exiguobacteraceae]; g__; s__" 
taxonomy_jamie_dict["b81ca62dc1c0e44c6f515c37a227d063"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Escherichia; s__coli" 
taxonomy_jamie_dict["c363a0062f4a061d01b9c52da4011aa1"]
  #"k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Escherichia; s__coli" 

#selbal like plot for coda-lasso
intervention_ba_2020.coda_pos <- intervention_ba_2020.results_codalasso$posCoefSelect
intervention_ba_2020.coda_neg <- intervention_ba_2020.results_codalasso$negCoefSelect
selbal_like_plot(pos.names = names(intervention_ba_2020.coda_pos),
                 neg.names = names(intervention_ba_2020.coda_neg),
                 Y = ba_intervention, X = table_intervention_ba_2020_1)

##Upset plot##
intervention_ba_2020.select <- list(selbal = intervention_ba.results_selbal$varSelect,
                                    clr_lasso = intervention_ba_2020.results_clrlasso$varSelect,
                                    coda_lasso = intervention_ba_2020.results_codalasso$varSelect)
intervention_ba_2020.select.upsetR <- fromList(intervention_ba_2020.select)
upset(as.data.frame(intervention_ba_2020.select.upsetR), main.bar.color = "gray36",
      sets.bar.color = color[c(1:2,5)], matrix.color = "gray36",
      order.by = "freq", empty.intersections = "on",
      queries = list(list(query = intersects, params = list("selbal"),
                          color = color[5], active = T),
                     list(query = intersects, params = list("clr_lasso"),
                          color = color[2], active = T),
                     list(query = intersects, params = list("coda_lasso"),
                          color = color[1], active = T)),
      text.scale = c(2, 2, 2, 2, 2, 2))


