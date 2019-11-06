library(ggplot2)
library(dplyr)

#df <- read.csv("rule_matches_vs_rule_non_matches_plotdata.csv", header = T)

rules <- read.csv("rule_matches_vs_rule_non_matches_plotdata.csv", header = T)

ALLrules <- rules[which(rules$comparison == "all rules"),]
LRrules <- rules[which(rules$comparison == "LR rules"),]
GLMrules <- LRrules[which(LRrules$algorithm == "GLM"),]
ARMrules <- rules[which(rules$comparison == "ARM rules"),]
ADASYNrules <- LRrules[which(LRrules$algorithm == "ADASYN"),]

mcr1rules <- read.csv("mcr_vs_nonmcr_self_vs_all_plotdata.csv", header = T)
clin1rules <- read.csv("clinical_vs_environmental_plotdata.csv", header = T)
clin2rules <- read.csv("clinical_vs_environmental_MCR_plotdata.csv", header = T)
chinrules <- read.csv("chinese_vs_non_chinese_plotdata.csv", header = T)

all_isolates <- mcr1rules[which(mcr1rules$set == "all"),]
non_mcr <- mcr1rules[which(mcr1rules$set == "non-mcr"),]
selfMCR <- mcr1rules[which(mcr1rules$set == "self"),]
self_no_mcr <- mcr1rules[which(mcr1rules$set == "self_without_mcr"),]
mcr_isolates <- mcr1rules[which(mcr1rules$set == "mcr"),]


#====================================
#            ALL RULES
#====================================


rule_match_all <- ALLrules[which(ALLrules$match == TRUE),]
lhs_match_all <- ALLrules[which(ALLrules$match_lhs == TRUE),]
non_lhs_match_all <- ALLrules[which(ALLrules$match_lhs == FALSE),]
non_match_all <- ALLrules[which(ALLrules$match == FALSE),]

all_df = rbind(rule_match_all,non_match_all)
match_df = rbind(rule_match_all,lhs_match_all)
non_df = rbind(non_match_all,non_lhs_match_all)

#====================================
#            ARM RULES
#====================================


rule_match_arm <- ARMrules[which(ARMrules$match == TRUE),]
lhs_match_arm <- ARMrules[which(ARMrules$match_lhs == TRUE),]
non_lhs_match_arm <- ARMrules[which(ARMrules$match_lhs == FALSE),]
non_match_arm <- ARMrules[which(ARMrules$match == FALSE),]

arm_df = rbind(rule_match_arm,non_match_arm)
matcharm_df = rbind(rule_match_arm,lhs_match_arm)
nonarm_df = rbind(non_match_arm,non_lhs_match_arm)

#====================================
#            GLM RULES
#====================================


rule_match_GLM <- GLMrules[which(GLMrules$match == TRUE),]
lhs_match_GLM <- GLMrules[which(GLMrules$match_lhs == TRUE),]
non_lhs_match_GLM <- GLMrules[which(GLMrules$match_lhs == FALSE),]
non_match_GLM <- GLMrules[which(GLMrules$match == FALSE),]

GLM_df = rbind(rule_match_GLM,non_match_GLM)
matchGLM_df = rbind(rule_match_GLM,lhs_match_GLM)
nonGLM_df = rbind(non_match_GLM,non_lhs_match_GLM)

#====================================
#            ADASYN RULES
#====================================


rule_match_ADASYN <- ADASYNrules[which(ADASYNrules$match == TRUE),]
lhs_match_ADASYN <- ADASYNrules[which(ADASYNrules$match_lhs == TRUE),]
non_lhs_match_ADASYN <- ADASYNrules[which(ADASYNrules$match_lhs == FALSE),]
non_match_ADASYN <- ADASYNrules[which(ADASYNrules$match == FALSE),]

ADASYN_df = rbind(rule_match_ADASYN,non_match_ADASYN)
matchADASYN_df = rbind(rule_match_ADASYN,lhs_match_ADASYN)
nonADASYN_df = rbind(non_match_ADASYN,non_lhs_match_ADASYN)

# lhs_match <- df[df$set == "rule lhs matches",]
# lhs_no_match <- df[df$set == "rule lhs non-matches",]
# rule_match <- df[df$set == "rule matches",]
# no_rule_match <- df[df$set == "rule non-matches",]

#lm <- ggplot(df %>% filter(set %in% c("rule lhs matches", "rule matches")),aes(x=pctself, color = set)) + geom_density(alpha = 0.2) + theme_classic() + xlab("Percentage of Individualistic Genes")

#====================================
#            PLOTTING
#====================================

chinese <- ggplot(chinrules,aes(x=(pctself*100),color = set, fill = set)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#E8BE4B", "#7E75A4"), 
                                                                                                                                                                       name="Location Comparison",
                                                                                                                                                                       breaks=c("chineseisolates", "nonchineseisolates"),
                                                                                                                                                                       labels=c("Isolates from China", "Isolates not from China")) + scale_fill_manual(values=c("#E8BE4B", "#7E75A4"), 
                                                                                                                                                                                                                                                        name="Location Comparison",
                                                                                                                                                                                                                                                        breaks=c("chineseisolates", "nonchineseisolates"),
                                                                                                                                                                                                                                                        labels=c("Isolates from China", "Isolates not from China"))
chinese + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

mcr1 <- ggplot(mcr1rules %>% filter(set %in% c("mcr", "self_without_mcr")),aes(x=(pctself*100), color = set, fill = set)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#CF4433", "#BE6EB7"), 
                                                                                                                                                                                                                      name="MCR Comparison",
                                                                                                                                                                                                                      breaks=c("mcr", "self_without_mcr"),
                                                                                                                                                                                                                      labels=c("mcr", "indiv wout mcr")) + scale_fill_manual(values=c("#CF4433", "#BE6EB7"), 
                                                                                                                                                                                                                                                                              name="MCR Comparison",
                                                                                                                                                                                                                                                                              breaks=c("mcr", "self_without_mcr"),
                                                                                                                                                                                                                                                                              labels=c("mcr", "indiv wout mcr"))
mcr1 + theme(legend.position = c(0.23, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

mcr2 <- ggplot(mcr1rules %>% filter(set %in% c("mcr", "non-mcr")),aes(x=(pctself*100), color = set, fill = set)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#CF4433", "#5EB9AD"), 
                                                                                                                                                                                                             name="MCR Comparison",
                                                                                                                                                                                                             breaks=c("mcr", "non-mcr"),
                                                                                                                                                                                                             labels=c("mcr", "no-mcr")) + scale_fill_manual(values=c("#CF4433", "#5EB9AD"), 
                                                                                                                                                                                                                                                             name="MCR Comparison",
                                                                                                                                                                                                                                                             breaks=c("mcr", "non-mcr"),
                                                                                                                                                                                                                                                             labels=c("mcr", "no-mcr"))
mcr2 + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

clinenv <- ggplot(clin1rules,aes(x=(pctself*100),color = set, fill = set)) + geom_density(alpha = 0.1) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#92B7D7", "#597CB5"), 
                                                                                                                                                                        name="Source Comparison",
                                                                                                                                                                        breaks=c("clinicalisolates", "environmentalisolates"),
                                                                                                                                                                        labels=c("Clinical", "Environmental")) + scale_fill_manual(values=c("#92B7D7", "#597CB5"), 
                                                                                                                                                                                                                                      name="Source Comparison",
                                                                                                                                                                                                                                      breaks=c("clinicalisolates", "environmentalisolates"),
                                                                                                                                                                                                                                      labels=c("Clinical", "Environmental"))
clinenv + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

clinenvmcr <- ggplot(clin2rules,aes(x=(pctself*100),color = set, fill = set)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")+ scale_color_manual(values=c("#D38472", "#CF4433"), 
                                                                                                                                                                          name="Source Comparison w/ MCR",
                                                                                                                                                                          breaks=c("clinicalmcrisolates", "environmentalmcrisolates"),
                                                                                                                                                                          labels=c("Clinical MCR", "Environmental MCR")) + scale_fill_manual(values=c("#D38472", "#CF4433"), 
                                                                                                                                                                                                                                             name="Source Comparison w/ MCR",
                                                                                                                                                                                                                                             breaks=c("clinicalmcrisolates", "environmentalmcrisolates"),
                                                                                                                                                                                                                                             labels=c("Clinical MCR", "Environmental MCR"))
clinenvmcr + theme(legend.position = c(0.27, 0.8), legend.title=element_text(size=7),legend.text=element_text(size=7))

#============================================

#               5 Generic Plots

#============================================

allDensity <- ggplot(ALLrules,aes(x=(pctself*100), xintercept = (mean(pctself)*100))) + geom_density(alpha = 1.0, fill = '#597CB5') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + geom_vline(aes(xintercept=mean(ALLrulez2$pctself)*100), linetype='dashed', color='#597CB5')
allDensity

non_mcrDensity <- ggplot(non_mcr,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#5EB9AD') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
non_mcrDensity

selfDensity <- ggplot(selfMCR,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#EA963C') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
selfDensity

self_noMCRDensity <- ggplot(self_no_mcr,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#BE6EB7') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
self_noMCRDensity

mcrDensity <- ggplot(mcr_isolates,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#CF4433') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
mcrDensity

#============================================

#         Individual Analysis Comparison

#============================================

#====================================
#            ARM RULES
#====================================

LHSrule_matchARMDensity <- ggplot(lhs_match_arm,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#9F1636') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
LHSrule_matchARMDensity

nonLHSrule_matchARMDensity <- ggplot(non_lhs_match_arm,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#C94245') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
nonLHSrule_matchARMDensity

rule_matchARMDensity <- ggplot(rule_match_arm,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#E26F5B') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
rule_matchARMDensity

nonrule_matchARMDensity <- ggplot(non_match_arm,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#ED9D89') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
nonrule_matchARMDensity


#====================================
#            GLM RULES
#====================================

LHSrule_matchGLMDensity <- ggplot(lhs_match_GLM,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#326A3C') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
LHSrule_matchGLMDensity

nonLHSrule_matchGLMDensity <- ggplot(non_lhs_match_GLM,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#508753') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
nonLHSrule_matchGLMDensity

rule_matchGLMDensity <- ggplot(rule_match_GLM,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#76A660') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
rule_matchGLMDensity

nonrule_matchGLMDensity <- ggplot(non_match_GLM,aes(x=(pctself*100))) + geom_density(alpha = 1.0, fill = '#A8BB67') + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")
nonrule_matchGLMDensity

#====================================
#            ALL RULES
#====================================
test <- non_df[is.na(non_df$match),]
test$match <- "missing"
non_df = rbind(non_match_all,test)

test <- match_df[is.na(match_df$match),]
test$match <- "missing"
match_df = rbind(rule_match_all,test)

all_match_vs_non  <- ggplot(all_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")+ scale_color_manual(values=c("#B68CAA", "#D3ADC9"), 
                                                                                                                                                                             name="ALL Matches vs. Non-Matches",
                                                                                                                                                                             breaks=c("TRUE", "FALSE"),
                                                                                                                                                                             labels=c("Rule Matches", "Rule Non-Matches"))+ scale_fill_manual(values=c("#B68CAA", "#D3ADC9"), 
                                                                                                                                                                                                                                               name="ALL Matches vs. Non-Matches",
                                                                                                                                                                                                                                               breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                                                                               labels=c("Rule Matches", "Rule Non-Matches")) 
all_match_vs_non + theme(legend.position = c(0.23, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

all_match_vs_lhs_match  <- ggplot(match_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#B68CAA", "#755076"), 
                                                                                                                                                                                       name="ALL Matches vs. LHS Matches",
                                                                                                                                                                                       breaks=c("TRUE", "missing"),
                                                                                                                                                                                       labels=c("Rule Matches", "LHS Matches")) + scale_fill_manual(values=c("#B68CAA", "#755076"), 
                                                                                                                                                                                                                                                     name="ALL Matches vs. LHS Matches",
                                                                                                                                                                                                                                                     breaks=c("TRUE", "missing"),
                                                                                                                                                                                                                                                     labels=c("Rule Matches", "LHS Matches"))
all_match_vs_lhs_match + theme(legend.position = c(0.23, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

all_non_vs_non_lhs  <- ggplot(non_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#D3ADC9", "#936D93"), 
                                                                                                                                                                               name="ALL Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                               breaks=c("FALSE", "missing"),
                                                                                                                                                                               labels=c("Rule Non-Matches", "LHS Non-Matches")) + scale_fill_manual(values=c("#D3ADC9", "#936D93"), 
                                                                                                                                                                                                                                                     name="ALL Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                                                                                     breaks=c("FALSE", "missing"),
                                                                                                                                                                                                                                                     labels=c("Rule Non-Matches", "LHS Non-Matches"))
all_non_vs_non_lhs + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

#====================================
#            ARM RULES
#====================================

test <- nonarm_df[is.na(nonarm_df$match),]
test$match <- "missing"
nonarm_df = rbind(non_match_arm,test)

test <- matcharm_df[is.na(matcharm_df$match),]
test$match <- "missing"
matcharm_df = rbind(rule_match_arm,test)

ARM_match_vs_non  <- ggplot(arm_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#E26F5B", "#ED9D89"), 
                                                                                                                                                                                     name="ARM Matches vs. Non-Matches",
                                                                                                                                                                                     breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                     labels=c("Rule Matches", "Rule Non-Matches"))+ scale_fill_manual(values=c("#E26F5B", "#ED9D89"), 
                                                                                                                                                                                                                                                       name="ARM Matches vs. Non-Matches",
                                                                                                                                                                                                                                                       breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                                                                                       labels=c("Rule Matches", "Rule Non-Matches"))
ARM_match_vs_non + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

ARM_match_vs_lhs_match  <- ggplot(matcharm_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#9F1636", "#E26F5B"), 
                                                                                                                                                                                                name="ARM Matches vs. LHS Matches",
                                                                                                                                                                                                breaks=c("TRUE", "missing"),
                                                                                                                                                                                                labels=c("Rule Matches", "LHS Matches"))+ scale_fill_manual(values=c("#E26F5B", "#9F1636"), 
                                                                                                                                                                                                                                                             name="ARM Matches vs. LHS Matches",
                                                                                                                                                                                                                                                             breaks=c("TRUE", "missing"),
                                                                                                                                                                                                                                                             labels=c("Rule Matches", "LHS Matches"))
ARM_match_vs_lhs_match + theme(legend.position = c(0.27, 0.8), legend.title=element_text(size=7),legend.text=element_text(size=7))

ARM_non_vs_non_lhs  <- ggplot(nonarm_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#ED9D89", "#C94245"), 
                                                                                                                                                                                           name="ARM Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                           breaks=c("FALSE", "missing"),
                                                                                                                                                                                           labels=c("Rule Non-Matches", "LHS Non-Matches"))+ scale_fill_manual(values=c("#ED9D89", "#C94245"), 
                                                                                                                                                                                                                                                                name="ARM Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                                                                                                breaks=c("FALSE", "missing"),
                                                                                                                                                                                                                                                                labels=c("Rule Non-Matches", "LHS Non-Matches"))
ARM_non_vs_non_lhs + theme(legend.position = c(0.43, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

#====================================
#            GLM RULES
#====================================

test <- nonGLM_df[is.na(nonGLM_df$match),]
test$match <- "missing"
nonGLM_df = rbind(non_match_GLM,test)

test <- matchGLM_df[is.na(matchGLM_df$match),]
test$match <- "missing"
matchGLM_df = rbind(rule_match_GLM,test)

GLM_match_vs_non  <- ggplot(GLM_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density")+ scale_color_manual(values=c("#76A660", "#A8BB67"), 
                                                                                                                                                                               name="GLM Matches vs. Non-Matches",
                                                                                                                                                                               breaks=c("TRUE", "FALSE"),
                                                                                                                                                                               labels=c("Rule Matches", "Rule Non-Matches")) + scale_fill_manual(values=c("#76A660", "#A8BB67"), 
                                                                                                                                                                                                                                                  name="GLM Matches vs. Non-Matches",
                                                                                                                                                                                                                                                  breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                                                                                  labels=c("Rule Matches", "Rule Non-Matches"))
GLM_match_vs_non + theme(legend.position = c(0.29, 0.8), legend.title=element_text(size=7),legend.text=element_text(size=7))

GLM_match_vs_lhs_match  <- ggplot(matchGLM_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#76A660", "#326A3C"), 
                                                                                                                                                                                          name="GLM Matches vs. LHS Matches",
                                                                                                                                                                                          breaks=c("TRUE", "missing"),
                                                                                                                                                                                          labels=c("Rule Matches", "LHS Matches")) + scale_fill_manual(values=c("#76A660", "#326A3C"), 
                                                                                                                                                                                                                                                        name="GLM Matches vs. LHS Matches",
                                                                                                                                                                                                                                                        breaks=c("TRUE", "missing"),
                                                                                                                                                                                                                                                        labels=c("Rule Matches", "LHS Matches"))
GLM_match_vs_lhs_match + theme(legend.position = c(0.29, 0.8), legend.title=element_text(size=7),legend.text=element_text(size=7))

GLM_non_vs_non_lhs  <- ggplot(nonGLM_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#A8BB67", "#508753"), 
                                                                                                                                                                                     name="GLM Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                     breaks=c("FALSE", "missing"),
                                                                                                                                                                                     labels=c("Rule Non-Matches", "LHS Non-Matches"))+ scale_fill_manual(values=c("#A8BB67", "#508753"), 
                                                                                                                                                                                                                                                          name="GLM Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                                                                                          breaks=c("FALSE", "missing"),
                                                                                                                                                                                                                                                          labels=c("Rule Non-Matches", "LHS Non-Matches"))
GLM_non_vs_non_lhs + theme(legend.position = c(0.43, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

#====================================
#            ADASYN RULES
#====================================

test <- nonADASYN_df[is.na(nonADASYN_df$match),]
test$match <- "missing"
nonADASYN_df = rbind(non_match_ADASYN,test)

test <- matchADASYN_df[is.na(matchADASYN_df$match),]
test$match <- "missing"
matchADASYN_df = rbind(rule_match_ADASYN,test)

ADASYN_match_vs_non  <- ggplot(ADASYN_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#749ABF", "#9BBFDC"), 
                                                                                                                                                                                                      name="ADASYN Matches vs. Non-Matches",
                                                                                                                                                                                                      breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                                      labels=c("Rule Matches", "Rule Non-Matches"))+ scale_fill_manual(values=c("#749ABF", "#9BBFDC"), 
                                                                                                                                                                                                                                                                        name="ADASYN Matches vs. Non-Matches",
                                                                                                                                                                                                                                                                        breaks=c("TRUE", "FALSE"),
                                                                                                                                                                                                                                                                        labels=c("Rule Matches", "Rule Non-Matches"))
ADASYN_match_vs_non + theme(legend.position = c(0.23, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

ADASYN_match_vs_lhs_match  <- ggplot(matchADASYN_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#749ABF", "#35577F"), 
                                                                                                                                                                                       name="ADASYN Matches vs. LHS Matches",
                                                                                                                                                                                       breaks=c("TRUE", "missing"),
                                                                                                                                                                                       labels=c("Rule Matches", "LHS Matches"))+ scale_fill_manual(values=c("#749ABF", "#35577F"), 
                                                                                                                                                                                                                                                    name="ADASYN Matches vs. LHS Matches",
                                                                                                                                                                                                                                                    breaks=c("TRUE", "missing"),
                                                                                                                                                                                                                                                    labels=c("Rule Matches", "LHS Matches"))
ADASYN_match_vs_lhs_match + theme(legend.position = c(0.23, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

ADASYN_non_vs_non_lhs  <- ggplot(nonADASYN_df,aes(x=(pctself*100),color = match, fill = match)) + geom_density(alpha = 0.2) + theme_classic() + xlim(0,100) + xlab("Percentage of Individualistic Genes") + ylab("Density") + scale_color_manual(values=c("#9BBFDC", "#5378A2"), 
                                                                                                                                                                                  name="ADASYN Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                  breaks=c("FALSE", "missing"),
                                                                                                                                                                                  labels=c("Rule Non-Matches", "LHS Non-Matches"))+ scale_fill_manual(values=c("#9BBFDC", "#5378A2"), 
                                                                                                                                                                                                                                                       name="ADASYN Non-Matches vs. LHS Non-Matches",
                                                                                                                                                                                                                                                       breaks=c("FALSE", "missing"),
                                                                                                                                                                                                                                                       labels=c("Rule Non-Matches", "LHS Non-Matches"))
ADASYN_non_vs_non_lhs + theme(legend.position = c(0.33, 0.8), legend.title=element_text(size=8),legend.text=element_text(size=8))

