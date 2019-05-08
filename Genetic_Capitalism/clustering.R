###############################################################
##    Cluster Analysis of AMR Genotype  Gains and Losses     ##
## By: Colby T. Ford, Ph.D., Gabriel Zenarosa, Ph.D.,        ##
##     Kevin Smith, John Williams, and Daniel Janies, Ph.D.  ##
###############################################################

library(dplyr)

## Get Genotype Counts
pdt <- unique(names(read.csv("e.coli.ids")))
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

## Read in Gain and Loss Data
data <- read.csv("gainloss_counts.csv")
data <- merge(data, counts, by = "Genotype")

# ## Prepare PWTs
# pwt <- data %>% filter(Genotype %in% c("blaEC", "blaEC-5"))
# pwt$rate <- pwt$Isolate_Count/(pwt$Gain-pwt$Loss)
# pwt$net_gain <- scale(pwt$Gain * pwt$rate)[,1]
# pwt$net_loss <- scale(pwt$Loss * pwt$rate)[,1]
# pwt$cluster <- NA
# pwt$category <- "PWT"


# ## Remove PWTs: blaEC & blaEC-5 (outliers)
# data <- data %>% filter(!Genotype %in% c("blaEC", "blaEC-5"))

## Calculate rate columns, net_gain, and net_loss
# data$pctdiff <- abs(data$Gain-data$Loss)/((data$Gain+data$Loss)/2)
# data$gain_rate <- data$Gain/data$Isolate_Count
# data$loss_rate <- data$Loss/data$Isolate_Count
data$rate <- data$Isolate_Count/(data$Gain-data$Loss)
data$net_gain <- (data$Gain * data$rate)
data$net_loss <- (data$Loss * data$rate)

############
## k-Means Analysis
library(mclust)

set.seed(1337)
meanclust <- Mclust(data[,6:7])

data$cluster <- factor(meanclust$classification)

cluster_counts <- data %>% 
  group_by(cluster) %>% 
  tally()

data <- data %>%
  mutate(category = case_when(cluster == 1 ~ "TBD",
                              cluster %in% c(2,4) ~ "GC",
                              cluster %in% c(3) ~ "SS"))
  # mutate(category = case_when(cluster == 1 ~ "TBD",
  #                             cluster %in% c(2,3,5) ~ "GC",
  #                             cluster %in% c(4,6) ~ "SS")) %>% 
  # rbind(pwt)

write.csv(data, "clusters.csv")

###############
## Plot Genotype Clusters
library(ggplot2)
library(plotly)


p <- ggplot(data = data,
             aes(x = Gain,
                 y = Loss,
                 color = cluster,
                 label = Genotype,
                 text = paste("Rate: ", rate,
                              "\nNet Gain: ", net_gain,
                              "\nNet Loss: ", net_loss))) + 
  geom_point(alpha = (2/3)) + 
  geom_label() +
  xlim(0,2500) +
  ylim(0,800)

ggplotly(p)
