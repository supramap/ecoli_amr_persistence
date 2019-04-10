##############################
##    Script to Generate    ##
##    Wilcoxon Test by      ##
##        AMR Group         ##
## By: Colby T. Ford, Ph.D. ##
##############################

library(dplyr)
data <- read.csv("Downloads/AMR_FunctionalMechanisms.csv")

######################
## GC/SS by Type
GCSSbyType <- data %>% group_by(GC.SS, Type, .drop = FALSE) %>% tally()
gc <- GCSSbyType %>% filter(GC.SS == "GC")
ss <- GCSSbyType %>% filter(GC.SS == "SS")

GCSSbyType_wt <- wilcox.test(gc$n,
                             ss$n,
                             alternative = "two.sided",
                             paired = TRUE)

GCSSbyType_wt$p.adjust <- p.adjust(GCSSbyType_wt$p.value,
                                   method = "bonferroni",
                                   n = nrow(gc))

######################
## GC/SS by Resistance Mechanism
GCSSbyRM <- data %>% group_by(GC.SS, Resistance.Mechanism, .drop = FALSE) %>% tally()
gc <- GCSSbyRM %>% filter(GC.SS == "GC")
ss <- GCSSbyRM %>% filter(GC.SS == "SS")

GCSSbyRM_wt <- wilcox.test(gc$n,
                           ss$n,
                           alternative = "two.sided",
                           paired = TRUE)

GCSSbyRM_wt$p.adjust <- p.adjust(GCSSbyRM_wt$p.value,
                                   method = "bonferroni",
                                   n = nrow(gc))

######################
## GC/SS by Gene
GCSSbyGene <- data %>% group_by(GC.SS, AMR.Gene, .drop = FALSE) %>% tally()
gc <- GCSSbyGene %>% filter(GC.SS == "GC")
ss <- GCSSbyGene %>% filter(GC.SS == "SS")

GCSSbyGene_wt <- wilcox.test(gc$n,
                             ss$n,
                             alternative = "two.sided",
                             paired = TRUE)

GCSSbyGene_wt$p.adjust <- p.adjust(GCSSbyGene_wt$p.value,
                                 method = "bonferroni",
                                 n = nrow(gc))
