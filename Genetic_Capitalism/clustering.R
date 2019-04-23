##############################
##    Cluster Analysis by   ##
##        AMR Genotype      ##
##      Gains and Losses    ##
## By: Colby T. Ford, Ph.D. ##
##############################
library(dplyr)

## Get Genotype Counts
pdt <- unique(names(read.csv("../MCR_Analysis/e.coli.ids")))
json.stream <- readRDS("TNTFileGenerator/e.coli.RDS")
id <- substr(json.stream[["ngout"]][["data"]][["content"]][["id"]], 19, 33)
fgenotypes <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(fgenotypes) <- id
rgenotypes <- fgenotypes
rgenotypes[!(names(fgenotypes) %in% pdt)] <- NULL

genotypesall <- unique(sort(unlist(rgenotypes)))
genobinv <- lapply(rgenotypes, function(x){(genotypesall %in% x) * 1})
gdf <- as.data.frame(matrix(unlist(genobinv), ncol = length(genotypesall), byrow = T))
names(gdf) <- genotypesall

counts <- as.data.frame(colSums(gdf))
colnames(counts) <- "Isolate_Count"
counts$Genotype <- rownames(counts) 


## Read in Data
data <- read.csv("gainloss_counts.csv")
#counts <- read.csv("IsolateCounts.csv")
data <- merge(data, counts, by = "Genotype")

## Remove blaEC & blaEC-5
data <- data %>% filter(!Genotype %in% c("blaEC", "blaEC-5"))

# data$pctdiff <- abs(data$Gain-data$Loss)/((data$Gain+data$Loss)/2)
# data$gain_rate <- data$Gain/data$Isolate_Count
# data$loss_rate <- data$Loss/data$Isolate_Count
data$rate <- data$Isolate_Count/(data$Gain-data$Loss)
data$net_gain <- scale(data$Gain * data$rate)
data$net_loss <- scale(data$Loss * data$rate)

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

cluster_counts <- data %>% 
  group_by(cluster) %>% 
  tally()

data <- data %>% mutate(cluster = case_when(cluster == 1 ~ "TBD",
                                            cluster == 2 ~ "GC",
                                            cluster == 3 ~ "GC",
                                            cluster == 4 ~ "SS",
                                            cluster == 5 ~ "GC",
                                            cluster == 6 ~ "SS"))
         
         
## Separate Clusters
# set.seed(1337)
# gainclust <- Mclust(data[,6])
# data$gaincluster <- factor(gainclust$classification)
# 
# lossclust <- Mclust(data[,7])
# data$losscluster <- factor(lossclust$classification)

# data <- data %>% 
#   mutate(gaincluster,
#          gaincluster = ifelse(gaincluster == 1,
#                               "Infrequently",
#                               "Frequently")) %>% 
#   mutate(losscluster,
#          losscluster = ifelse(losscluster == 1,
#                               "Infrequently",
#                               "Frequently"))
# 
# data$cluster <- paste0("G: ",
#                        data$gaincluster,
#                        " | L: "
#                        , data$losscluster)
# 
# ## Define GC, SS, and TBD Groups
# data$GC.SS <- "SS"
# data$GC.SS[data$cluster == "G: Frequently | L: Infrequently"] <- "GC"
# data$GC.SS[data$cluster == "G: Infrequently | L: Infrequently"] <- "TBD"
# 
# cluster_counts <- data %>% group_by(GC.SS) %>% tally()

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
# cluster_counts <- data %>% group_by(cluster) %>% tally()

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
             aes(x=Gain,
                 y=Loss,
                 color=cluster,
                 label=Genotype)) + 
  geom_point() + 
  geom_label() #+
  # ylim(0,1) +
  # xlim(0,1)
ggplotly(lp)


