##############################
##    Cluster Analysis by   ##
##        AMR Genotype      ##
##      Gains and Losses    ##
## By: Colby T. Ford, Ph.D. ##
##############################

## Read in Data
data <- read.csv("gainloss_counts.csv")
data$pctdiff <- abs(data$Gain-data$Loss)/((data$Gain+data$Loss)/2)

############
## k-Means Analysis
# set.seed(1337)
# meanclust <- kmeans(scale(data[,1:2]),
#                     centers = 5,
#                     iter.max = 10000)
# 
# data$cluster <- factor(meanclust$cluster)
# centers <- as.data.frame(meanclust$centers)

library(mclust)
meanclust <- Mclust(data[,1:2])
#meanclust <- Mclust(data[,3])

data$cluster <- factor(meanclust$classification)

counts <- data %>% group_by(cluster) %>% tally()

## Separate Clusters
gainclust <- Mclust(data[,1])
data$gaincluster <- factor(gainclust$classification)
lossclust <- Mclust(data[,2])
data$losscluster <- factor(lossclust$classification)

data$cluster <- paste0(data$gaincluster, "|", data$losscluster)

############
## k-Median Analysis
# library(Gmedian)
# set.seed(1337)
# medclust <- kGmedian(data[,1:2],
#                      ncenters = 5)
# data$cluster <- factor(medclust$cluster)
# centers <- as.data.frame(medclust$centers)
# 
# counts <- data %>% group_by(cluster) %>% tally()

###############
## Plot Genotype Clusters
library(ggplot2)
library(plotly)

p <- ggplot(data = data, aes(x=Gain, y=Loss, color=cluster, label=Genotype)) + 
  geom_point() + 
  geom_label() +
  geom_vline(xintercept = 20) +
  geom_hline(yintercept = 17) +
  ylim(0,2500)

p <- ggplot(data = data, aes(x=log(Gain), y=log(Loss), color=cluster, label=Genotype)) + 
  geom_point() + 
  geom_label()

ggplotly(p)


## Define GC, SS, and TBD Groups
data$GC.SS <- "SS"
data$GC.SS[data$cluster == "2|1"] <- "GC"
data$GC.SS[data$cluster == "1|1"] <- "TBD"

write.csv(data, "clusters.csv")
