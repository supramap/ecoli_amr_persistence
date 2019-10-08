############################################################
##        Script for Evaludating Distributions of         ##
##            of Individualistic Gene Content             ##
##                                                        ##
##        By: Kevin Smith, Colby T. Ford, Ph.D.,          ##
##         Gabriel Zenarosa, Ph.D., David Brown,          ##
##               and Daniel Janies, Ph.D.                 ##
############################################################

## Load in Libraries
library(dplyr)
library(ggplot2)
library(readr)

######################
## MCR vs. non-MCR vs. Individualistic
mcr_vs_nonmcr_self_vs_all_plotdata <- read.csv("mcr_vs_nonmcr_self_vs_all_plotdata.csv")

### MCR vs Non-MCR_Selfish (Individualistic (without MCR))
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "self_without_mcr")),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "self_without_mcr")),
            paired = FALSE, alternative = "greater")



### MCR vs Non-MCR
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "non-mcr")),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "non-mcr")),
            paired = FALSE, alternative = "greater")


########################

## Chinese vs. Non-Chinese Isolates
chinese_vs_non_chinese_plotdata <- read.csv("chinese_vs_non_chinese_plotdata.csv")

wilcox.test(pctself ~ set,
            data = chinese_vs_non_chinese_plotdata %>%
              filter(set %in% c("chineseisolates", "nonchineseisolates")),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = chinese_vs_non_chinese_plotdata %>%
              filter(set %in% c("chineseisolates", "nonchineseisolates")),
            paired = FALSE, alternative = "greater")

## Clinical vs. Environmental Isolates
clinical_vs_environmental_plotdata <- read.csv("clinical_vs_environmental_plotdata.csv")

wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_plotdata,
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_plotdata,
            paired = FALSE, alternative = "less")

## Clinical vs. Environmental MCR Isolates
clinical_vs_environmental_MCR_plotdata <- read.csv("clinical_vs_environmental_MCR_plotdata.csv")

wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_MCR_plotdata,
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_MCR_plotdata,
            paired = FALSE, alternative = "greater")

########################
## Rule vs. Isolate Set Comparison
rule_matches_vs_rule_non_matches_plotdata <- readRDS("rule_matches_vs_rule_non_matches_plotdata.RDS")

### SET COMPARISONS - STATS TESTS (ALL)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "All"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "greater")

## lhs Rule matches vs. non-matches (without MCR)
# wilcox.test(pctself ~ set,
#             data = rule_matches_vs_rule_non_matches_plotdata %>% 
#               filter(set %in% c("rule lhs matches", "rule lhs non-matches"),
#                      algorithm == "All"),
#             paired = FALSE)

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "All"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "All"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (ARM)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "apriori"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "apriori"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "apriori"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (GLM)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "GLM"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "GLM"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "GLM"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (ADASYN)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "ADASYN"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE)

wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "less")
