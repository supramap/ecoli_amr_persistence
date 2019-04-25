###############################################################
##    Distribution Tests by AMR Group & Resistance Mechanism ##
## By: Colby T. Ford, Ph.D., Gabriel Zenarosa, Ph.D.,        ##
##     Kevin Smith, John Williams, and Daniel Janies, Ph.D.  ##
###############################################################

library(dplyr)
library(ggplot2)
library(plotly)
library(ggsci)

## Load in resistance mechanism data
data <- read.csv("AMR_FunctionalMechanisms.csv")

## Remove blaEC & blaEC-5
data <- data %>% filter(!AMR.Gene %in% c("blaEC", "blaEC-5"))

## Read in new clusters
newgcss <- read.csv("clusters.csv")
data$GC.SS <- newgcss$cluster

######################
## GC/SS by Type
GCSSbyType <- data %>% group_by(GC.SS, Type, .drop=FALSE) %>% tally()

gc <- GCSSbyType %>% filter(GC.SS == "GC")
gc$pct <- gc$n/sum(gc$n)

ss <- GCSSbyType %>% filter(GC.SS == "SS")
ss$pct <- ss$n/sum(ss$n)

tbd <- GCSSbyType %>% filter(GC.SS == "TBD")
tbd$pct <- tbd$n/sum(tbd$n)


## Fisher Test

## For all 3 categories: GC, SS, TBD
GCSSbyType_ft_3 <- fisher.test(matrix(c(gc$n,
                                        ss$n,
                                        tbd$n),
                                      3,
                                      6,
                                      byrow=TRUE,
                                      dimnames = list(classification=c("GC",
                                                                       "SS",
                                                                       "TBD"),
                                                      mechanism=c("unknown",
                                                                  "alteration",
                                                                  "efflux",
                                                                  "inactiviation",
                                                                  "protection",
                                                                  "replacement"))),
                               simulate.p.value = TRUE)

## For 2 categories: GC, SS
GCSSbyType_ft_2 <- fisher.test(matrix(c(gc$n,
                                        ss$n),
                                      2,
                                      6,
                                      byrow=TRUE,
                                      dimnames = list(classification=c("GC",
                                                                       "SS"),
                                                      mechanism=c("unknown",
                                                                  "alteration",
                                                                  "efflux",
                                                                  "inactiviation",
                                                                  "protection",
                                                                  "replacement"))),
                               simulate.p.value = TRUE)


## Plotting
p <- ggplot(data, aes(Type,  fill = GC.SS)) +
  geom_bar(position="dodge") +
  scale_fill_npg()

ggplotly(p)
