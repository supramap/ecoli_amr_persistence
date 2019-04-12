##############################
##    Script to Generate    ##
##  Distribution Tests by   ##
##        AMR Group         ##
## By: Colby T. Ford, Ph.D. ##
##############################

library(dplyr)
library(ggplot2)
library(plotly)
library(ggsci)

data <- read.csv("AMR_FunctionalMechanisms.csv")
newgcss <- read.csv("clusters.csv")
data$GC.SS <- newgcss$GC.SS

######################
## GC/SS by Type
GCSSbyType <- data %>% group_by(GC.SS, Type, .drop = FALSE) %>% tally()

gc <- GCSSbyType %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)

ss <- GCSSbyType %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

# GCSSbyType_wt <- wilcox.test(gc$n,
#                              ss$n,
#                              alternative = "two.sided",
#                              paired = TRUE)

GCSSbyType_ks <- ks.test(gc$n,
                         ss$n,
                         alternative = "two.sided",
                         exact = TRUE)

p <- ggplot(data, aes(Type, fill = GC.SS)) +
  geom_bar(position="dodge") +
  scale_fill_npg()

ggplotly(p)

######################
## GC/SS by Resistance Mechanism
GCSSbyRM <- data %>% group_by(GC.SS, Resistance.Mechanism, .drop = FALSE) %>% tally()
gc <- GCSSbyRM %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)
ss <- GCSSbyRM %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

# GCSSbyRM_wt <- wilcox.test(gc$n,
#                            ss$n,
#                            alternative = "two.sided",
#                            paired = TRUE)

GCSSbyRM_ks <- ks.test(gc$pct,
                       ss$pct,
                       alternative = "two.sided",
                       exact = TRUE)

p <- ggplot(data, aes(Resistance.Mechanism, fill = GC.SS)) +
  geom_bar(position="dodge") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_npg() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

ggplotly(p)

######################
## GC/SS by Gene
GCSSbyGene <- data %>% group_by(GC.SS, AMR.Gene, .drop = FALSE) %>% tally()

gc <- GCSSbyGene %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)

ss <- GCSSbyGene %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

# GCSSbyGene_wt <- wilcox.test(gc$n,
#                              ss$n,
#                              alternative = "two.sided",
#                              paired = TRUE)

GCSSbyGene_ks <- ks.test(gc$n,
                         ss$n,
                         alternative = "two.sided",
                         paired = TRUE)
