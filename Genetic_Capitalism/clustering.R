##############################
##    Cluster Analysis by   ##
##        AMR Genotype      ##
##      Gains and Losses    ##
## By: Colby T. Ford, Ph.D. ##
##############################
library(dplyr)
## Read in Data
data <- read.csv("gainloss_counts.csv")
counts <- read.csv("IsolateCounts.csv")
data <- merge(data, counts, by = "Genotype")
data$pctdiff <- abs(data$Gain-data$Loss)/((data$Gain+data$Loss)/2)
data$gain_rate <- data$Gain/data$Isolate_Count
data$loss_rate <- data$Loss/data$Isolate_Count

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
## Combined Clusters
set.seed(1337)
meanclust <- Mclust(data[,6:7])
#meanclust <- Mclust(data[,3])

data$cluster <- factor(meanclust$classification)

counts <- data %>% group_by(cluster) %>% tally()

## Separate Clusters
set.seed(1337)
gainclust <- Mclust(data[,6])
data$gaincluster <- factor(gainclust$classification)

lossclust <- Mclust(data[,7])
data$losscluster <- factor(lossclust$classification)

data <- data %>% 
  mutate(gaincluster,
         gaincluster = ifelse(gaincluster == 1,
                              "Infrequently",
                              "Frequently")) %>% 
  mutate(losscluster,
         losscluster = ifelse(losscluster == 1,
                              "Infrequently",
                              "Frequently"))

data$cluster <- paste0("G: ",
                       data$gaincluster,
                       " | L: "
                       , data$losscluster)

## Define GC, SS, and TBD Groups
data$GC.SS <- "SS"
data$GC.SS[data$cluster == "G: Frequently | L: Infrequently"] <- "GC"
data$GC.SS[data$cluster == "G: Infrequently | L: Infrequently"] <- "TBD"

counts <- data %>% group_by(GC.SS) %>% tally()

write.csv(data, "clusters.csv")

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

p <- ggplot(data = data, aes(x=Gain,
                             y=Loss,
                             color=GC.SS,
                             label=Genotype)) + 
  geom_point() + 
  geom_label() +
  #geom_vline(xintercept = 20) +
  #geom_hline(yintercept = 17) +
  ylim(0,2500)

ggplotly(p)

lp <- ggplot(data = data,
            aes(x=log(Gain),
                y=log(Loss),
                color=GC.SS,
                label=Genotype)) + 
  geom_point() + 
  geom_label()

ggplotly(lp)


lp <- ggplot(data = data,
             aes(x=gain_rate,
                 y=loss_rate,
                 color=cluster,
                 label=Genotype)) + 
  geom_point() + 
  geom_label() #+
  # ylim(0,1) +
  # xlim(0,1)
ggplotly(lp)


