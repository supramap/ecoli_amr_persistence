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

GCSSbyType_wt <- wilcox.test(gc$n, ss$n)


######################
## GC/SS by Resistance Mechanism
GCSSbyRM <- data %>% group_by(GC.SS, Resistance.Mechanism, .drop = FALSE) %>% tally()
gc <- GCSSbyRM %>% filter(GC.SS == "GC")
ss <- GCSSbyRM %>% filter(GC.SS == "SS")

GCSSbyRM_wt <- wilcox.test(gc$n, ss$n)
