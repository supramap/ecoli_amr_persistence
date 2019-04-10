##############################
##    Script to Generate    ##
##    Wilcoxon Test by      ##
##        AMR Group         ##
## By: Colby T. Ford, Ph.D. ##
##############################

library(dplyr)
data <- read.csv("AMR_FunctionalMechanisms.csv")

######################
## GC/SS by Type
GCSSbyType <- data %>% group_by(GC.SS, Type, .drop = FALSE) %>% tally()

gc <- GCSSbyType %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)

ss <- GCSSbyType %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

GCSSbyType_wt <- wilcox.test(gc$n,
                             ss$n,
                             alternative = "two.sided",
                             paired = TRUE)

GCSSbyType_ks <- ks.test(gc$n,
                         ss$n,
                         alternative = "two.sided",
                         exact = TRUE)

######################
## GC/SS by Resistance Mechanism
GCSSbyRM <- data %>% group_by(GC.SS, Resistance.Mechanism, .drop = FALSE) %>% tally()
gc <- GCSSbyRM %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)
ss <- GCSSbyRM %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

GCSSbyRM_wt <- wilcox.test(gc$n,
                           ss$n,
                           alternative = "two.sided",
                           paired = TRUE)

GCSSbyRM_ks <- ks.test(gc$pct,
                       ss$pct,
                       alternative = "two.sided",
                       exact = TRUE)


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
