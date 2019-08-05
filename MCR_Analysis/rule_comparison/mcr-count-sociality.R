library(tidyr)
library(dplyr)

######################
## For Association Rules
rules <- readRDS("../mcr_rules_all.RDS")
gene <- read.csv("../AMR_FunctionalMechanisms.csv", header=T)

#Identify genes are their respective socialities
for (i in 1:nrow(gene)) {
  coop = ifelse(gene$Sociality == "Cooperative", paste(gene$AMR.Gene), NA)
  self = ifelse(gene$Sociality == "Selfish", paste(gene$AMR.Gene), NA)
  gene$coop <- coop
  gene$self <- self
}

#Makes vectors contaning genes in their categories
coop_vector <- c(gene$coop[!is.na(gene$coop)])
self_vector <- c(gene$self[!is.na(gene$self)])

#Create new data frame for which to individually count genes on socaility metric
newruleslhs <- rules %>% tidyr::separate(lhs,sep=",",c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"),extra='drop')

#Count Genes and identify which genes are not counted to verify sociality metrics
for (i in 1:nrow(newruleslhs)) {
  newruleslhs$"1gene" <- ifelse(newruleslhs$"1" %in% coop_vector,1,0) + ifelse(newruleslhs$"1" %in% self_vector,1,0)
  newruleslhs$"2gene" <- ifelse(newruleslhs$"2" %in% coop_vector,1,0) + ifelse(newruleslhs$"2" %in% self_vector,1,0)
  newruleslhs$"3gene" <- ifelse(newruleslhs$"3" %in% coop_vector,1,0) + ifelse(newruleslhs$"3" %in% self_vector,1,0)
  newruleslhs$"4gene" <- ifelse(newruleslhs$"4" %in% coop_vector,1,0) + ifelse(newruleslhs$"4" %in% self_vector,1,0)
  newruleslhs$"5gene" <- ifelse(newruleslhs$"5" %in% coop_vector,1,0) + ifelse(newruleslhs$"5" %in% self_vector,1,0)
  newruleslhs$"6gene" <- ifelse(newruleslhs$"6" %in% coop_vector,1,0) + ifelse(newruleslhs$"6" %in% self_vector,1,0)
  newruleslhs$"7gene" <- ifelse(newruleslhs$"7" %in% coop_vector,1,0) + ifelse(newruleslhs$"7" %in% self_vector,1,0)
  newruleslhs$"8gene" <- ifelse(newruleslhs$"8" %in% coop_vector,1,0) + ifelse(newruleslhs$"8" %in% self_vector,1,0)
  newruleslhs$"9gene" <- ifelse(newruleslhs$"9" %in% coop_vector,1,0) + ifelse(newruleslhs$"9" %in% self_vector,1,0)
  newruleslhs$"10gene" <- ifelse(newruleslhs$"10" %in% coop_vector,1,0) + ifelse(newruleslhs$"10" %in% self_vector,1,0)
  newruleslhs$"11gene" <- ifelse(newruleslhs$"11" %in% coop_vector,1,0) + ifelse(newruleslhs$"11" %in% self_vector,1,0)
  newruleslhs$"12gene" <- ifelse(newruleslhs$"12" %in% coop_vector,1,0) + ifelse(newruleslhs$"12" %in% self_vector,1,0)
  newruleslhs$"13gene" <- ifelse(newruleslhs$"13" %in% coop_vector,1,0) + ifelse(newruleslhs$"13" %in% self_vector,1,0)
  newruleslhs$"14gene" <- ifelse(newruleslhs$"14" %in% coop_vector,1,0) + ifelse(newruleslhs$"14" %in% self_vector,1,0)
  newruleslhs$"15gene" <- ifelse(newruleslhs$"15" %in% coop_vector,1,0) + ifelse(newruleslhs$"15" %in% self_vector,1,0)
  newruleslhs$"16gene" <- ifelse(newruleslhs$"16" %in% coop_vector,1,0) + ifelse(newruleslhs$"16" %in% self_vector,1,0)
  newruleslhs$"17gene" <- ifelse(newruleslhs$"17" %in% coop_vector,1,0) + ifelse(newruleslhs$"17" %in% self_vector,1,0)
  
  numcoop <- ifelse(newruleslhs$"1" %in% coop_vector,1,0) + ifelse(newruleslhs$"2" %in% coop_vector,1,0)+ ifelse(newruleslhs$"3" %in% coop_vector,1,0)+ ifelse(newruleslhs$"4" %in% coop_vector,1,0)+ ifelse(newruleslhs$"5" %in% coop_vector,1,0)+ ifelse(newruleslhs$"6" %in% coop_vector,1,0)+ ifelse(newruleslhs$"7" %in% coop_vector,1,0)+ ifelse(newruleslhs$"8" %in% coop_vector,1,0)+ ifelse(newruleslhs$"9" %in% coop_vector,1,0)+ ifelse(newruleslhs$"10" %in% coop_vector,1,0)+ ifelse(newruleslhs$"11" %in% coop_vector,1,0)+ ifelse(newruleslhs$"12" %in% coop_vector,1,0)+ ifelse(newruleslhs$"13" %in% coop_vector,1,0)+ ifelse(newruleslhs$"14" %in% coop_vector,1,0)+ ifelse(newruleslhs$"15" %in% coop_vector,1,0)+ ifelse(newruleslhs$"16" %in% coop_vector,1,0)+ ifelse(newruleslhs$"17" %in% coop_vector,1,0)
  newruleslhs$numcoop <- numcoop
  rules$numcoop <- numcoop
  
  numself <- ifelse(newruleslhs$"1" %in% self_vector,1,0) + ifelse(newruleslhs$"2" %in% self_vector,1,0)+ ifelse(newruleslhs$"3" %in% self_vector,1,0)+ ifelse(newruleslhs$"4" %in% self_vector,1,0)+ ifelse(newruleslhs$"5" %in% self_vector,1,0)+ ifelse(newruleslhs$"6" %in% self_vector,1,0)+ ifelse(newruleslhs$"7" %in% self_vector,1,0)+ ifelse(newruleslhs$"8" %in% self_vector,1,0)+ ifelse(newruleslhs$"9" %in% self_vector,1,0)+ ifelse(newruleslhs$"10" %in% self_vector,1,0)+ ifelse(newruleslhs$"11" %in% self_vector,1,0)+ ifelse(newruleslhs$"12" %in% self_vector,1,0)+ ifelse(newruleslhs$"13" %in% self_vector,1,0)+ ifelse(newruleslhs$"14" %in% self_vector,1,0)+ ifelse(newruleslhs$"15" %in% self_vector,1,0)+ ifelse(newruleslhs$"16" %in% self_vector,1,0)+ ifelse(newruleslhs$"17" %in% self_vector,1,0)
  newruleslhs$numself <- numself
  rules$numself <- numself
}

#==========RMSD Calculations===========================================================
library(dplyr)
all_rules <- readRDS("../mcr_rules_all.RDS")
val_rules <- readRDS("../mcr_validation_rules.RDS")

#Join excluding column-list called "rule" (column 8) because of error and irrelevancy
inner_rules <- inner_join(all_rules[,-8], val_rules[,-8], by="lhs")
left_rules <- left_join(all_rules[,-8], val_rules[,-8], by="lhs")

#Set non-matching sets conf,supp,lift to 0
left_rules$support.y[is.na(left_rules$support.y)] <- 0
left_rules$confidence.y[is.na(left_rules$confidence.y)] <- 1
# left_rules$lift.y[is.na(left_rules$lift.y)] <- 1

#RMSD Calculations
mse_inner_support <- mean((inner_rules$support.x - inner_rules$support.y)^2)
mse_left_support <- mean((left_rules$support.x - left_rules$support.y)^2)
rmsd_inner_support <- sqrt((sum((abs(inner_rules$support.x - inner_rules$support.y))^2))/nrow(inner_rules))
rmsd_left_support <- sqrt((sum((abs(left_rules$support.x - left_rules$support.y))^2))/nrow(left_rules))

mse_inner_confidence <- mean((inner_rules$confidence.x - inner_rules$confidence.y)^2)
mse_left_confidence <- mean((left_rules$confidence.x - left_rules$confidence.y)^2)
rmsd_inner_confidence <- sqrt((sum((abs(inner_rules$confidence.x - inner_rules$confidence.y))^2))/nrow(inner_rules))
rmsd_left_confidence <- sqrt((sum((abs(left_rules$confidence.x - left_rules$confidence.y))^2))/nrow(left_rules))

# rmsd_inner_lift <- sqrt((sum((abs(inner_rules$lift.x - inner_rules$lift.y))^2))/nrow(inner_rules))
# rmsd_left_lift <- sqrt((sum((abs(left_rules$lift.x - left_rules$lift.y))^2))/nrow(left_rules))

###############################
## Coop/Self Ratio for MCR vs. Not
## Distribution of coopertivity of MCR-containing isolates
## Look in LR script to line 51

## Cols w/ coop-prefix/everything else (for mcr containing vs. not)
library(jsonlite)
library(stringr)
library(glmnet)
library(caret)
#library(arules)

# Run these *nix commands
#sed 's/'\''//g' PDG000000004.1024.reference_target.tree.newick | sed 's/:-*[0-9]\.*[0-9]*\(e-[0-9]*\)*//g' | sed 's/,'\('/'\('/g' | sed 's/,/ /g' > e.coli.paren
#sed 's/'\)'*'\('*'\ '*PDT/,PDT/g' e.coli.paren | sed 's/^,//g' | sed 's/'\)'*;$//g' > e.coli.ids
#---------
pdt <- unique(names(read.csv("../e.coli.ids")))

#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=1000000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Casc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=asc)")
#saveRDS(json.stream, file = "e.coli.RDS")
json.stream <- readRDS("../e.coli.RDS")

id <- substr(json.stream[["ngout"]][["data"]][["content"]][["id"]], 19, 33)
creation_date_time <- as.POSIXct(json.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")
collection_year <- as.numeric(substr(json.stream[["ngout"]][["data"]][["content"]][["collection_date"]], 1, 4))
collection_year[is.na(collection_year)] <- as.numeric(format(creation_date_time[is.na(collection_year)], "%Y"))
location <- json.stream[["ngout"]][["data"]][["content"]][["geo_loc_name"]]
isolation_type <- json.stream[["ngout"]][["data"]][["content"]][["epi_type"]]
isolation_source <- json.stream[["ngout"]][["data"]][["content"]][["isolation_source"]]
fdf <- data.frame(id, creation_date_time, collection_year, location, isolation_type, isolation_source)
rdf <- fdf[fdf$id %in% pdt,]
fgenotypes <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(fgenotypes) <- id
rgenotypes <- fgenotypes
rgenotypes[!(names(fgenotypes) %in% pdt)] <- NULL
#fgenotypes[fgenotypes == "NULL"] <- NULL
#rgenotypes[rgenotypes == "NULL"] <- NULL
genotypes <- fgenotypes[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]
genoid <- id[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]

# ---------
# Genotypes
# ---------
genotypesall <- unique(sort(unlist(genotypes)))
genobinv <- lapply(genotypes, function(x){(genotypesall %in% x) * 1})
gdf <- as.data.frame(matrix(unlist(genobinv), ncol = length(genotypesall), byrow = T))
names(gdf) <- genotypesall
gdf.mcr <- gdf[,-grep("mcr", names(gdf))]
gdf.mcr$mcr <- (rowSums(gdf[,grep("mcr", names(gdf))]) > 0) * 1
# gdf.mcr$mcr <- as.factor(gdf.mcr$mcr)

allisolates <- gdf.mcr
onlymcrisolates <- gdf.mcr %>% filter(mcr == 1)
onlyblaisolates <- gdf.mcr %>% filter_at(vars(starts_with("bla")), any_vars(. == 1))
onlyrmtisolates <- gdf.mcr %>% filter_at(vars(starts_with("rmt")), any_vars(. == 1))
onlyermisolates <- gdf.mcr %>% filter_at(vars(starts_with("erm")), any_vars(. == 1))


## Identify genes are their respective socialities
gene <- read.csv("../AMR_FunctionalMechanisms.csv", header=T)

for (i in 1:nrow(gene)) {
  coop = ifelse(gene$Sociality == "Cooperative", paste(gene$AMR.Gene), NA)
  self = ifelse(gene$Sociality == "Selfish", paste(gene$AMR.Gene), NA)
  unkn = ifelse(gene$Sociality == "", paste(gene$AMR.Gene), NA)
  gene$coop <- coop
  gene$self <- self
  gene$unkn <- unkn
}

## Makes vectors contaning genes in their categories
coop_vector <- c(gene$coop[!is.na(gene$coop)])
self_vector <- c(gene$self[!is.na(gene$self)])
unkn_vector <- c(gene$self[!is.na(gene$unkn)])

coop_vector_ors <- paste(coop_vector, collapse = "|")
self_vector_ors <- paste(self_vector, collapse = "|")
unkn_vector_ors <- paste(unkn_vector, collapse = "|")
all_vector_ors <- paste0(coop_vector_ors,"|",self_vector_ors,"|",unkn_vector_ors)

## Get isolates with at least 1 individualist gene
onlyselfisolates <- gdf.mcr %>% filter_at(vars(matches(self_vector_ors)), any_vars(. == 1))

## Look at ratios overall vs. only MCR vs only ...

allisolates <- allisolates %>% 
  mutate(set = "all",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(allisolates$pctself)


onlyselfisolates <- onlyselfisolates %>% 
  mutate(set = "self",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(onlyindividualisticisolates$pctself)

onlymcrisolates <- onlymcrisolates %>% 
  mutate(set = "mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(onlymcrisolates$pctself)


onlyblaisolates <- onlyblaisolates %>% 
  mutate(set = "bla",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(onlyblaisolates$pctself)


onlyrmtisolates <- onlyrmtisolates %>% 
  mutate(set = "rmt",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(onlyrmtisolates$pctself)

onlyermisolates <- onlyermisolates %>% 
  mutate(set = "erm",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count)

hist(onlyermisolates$pctself)

######################
## For Regression Models
rules <- read.csv("../mcr_lr_model_oddsratios.csv") %>%
  filter(OddsRatio_adasyn >= 1) %>% 
  mutate(GeneSet = stringr::str_replace_all(.$GeneSet,"`",""))

rules$genevector <- strsplit(rules$GeneSet,":")
rules$size <- lengths(rules$genevector)+1

gene <- read.csv("../AMR_FunctionalMechanisms.csv", header=T)

#Identify genes are their respective socialities
for (i in 1:nrow(gene)) {
  coop = ifelse(gene$Sociality == "Cooperative", paste(gene$AMR.Gene), NA)
  self = ifelse(gene$Sociality == "Selfish", paste(gene$AMR.Gene), NA)
  gene$coop <- coop
  gene$self <- self
}

#Makes vectors contaning genes in their categories
coop_vector <- c(gene$coop[!is.na(gene$coop)])
self_vector <- c(gene$self[!is.na(gene$self)])

#Create new data frame for which to individually count genes on socaility metric
newruleslhs <- rules %>%  
  tidyr::separate(GeneSet,sep=":",
                  c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),
                  extra='drop')
                                         

#Count Genes and identify which genes are not counted to verify sociality metrics
for (i in 1:nrow(newruleslhs)) {
  newruleslhs$"1gene" <- ifelse(newruleslhs$"1" %in% coop_vector,1,0) + ifelse(newruleslhs$"1" %in% self_vector,1,0)
  newruleslhs$"2gene" <- ifelse(newruleslhs$"2" %in% coop_vector,1,0) + ifelse(newruleslhs$"2" %in% self_vector,1,0)
  newruleslhs$"3gene" <- ifelse(newruleslhs$"3" %in% coop_vector,1,0) + ifelse(newruleslhs$"3" %in% self_vector,1,0)
  newruleslhs$"4gene" <- ifelse(newruleslhs$"4" %in% coop_vector,1,0) + ifelse(newruleslhs$"4" %in% self_vector,1,0)
  newruleslhs$"5gene" <- ifelse(newruleslhs$"5" %in% coop_vector,1,0) + ifelse(newruleslhs$"5" %in% self_vector,1,0)
  newruleslhs$"6gene" <- ifelse(newruleslhs$"6" %in% coop_vector,1,0) + ifelse(newruleslhs$"6" %in% self_vector,1,0)
  newruleslhs$"7gene" <- ifelse(newruleslhs$"7" %in% coop_vector,1,0) + ifelse(newruleslhs$"7" %in% self_vector,1,0)
  newruleslhs$"8gene" <- ifelse(newruleslhs$"8" %in% coop_vector,1,0) + ifelse(newruleslhs$"8" %in% self_vector,1,0)
  newruleslhs$"9gene" <- ifelse(newruleslhs$"9" %in% coop_vector,1,0) + ifelse(newruleslhs$"9" %in% self_vector,1,0)
  newruleslhs$"10gene" <- ifelse(newruleslhs$"10" %in% coop_vector,1,0) + ifelse(newruleslhs$"10" %in% self_vector,1,0)
  newruleslhs$"11gene" <- ifelse(newruleslhs$"11" %in% coop_vector,1,0) + ifelse(newruleslhs$"11" %in% self_vector,1,0)
  newruleslhs$"12gene" <- ifelse(newruleslhs$"12" %in% coop_vector,1,0) + ifelse(newruleslhs$"12" %in% self_vector,1,0)
  newruleslhs$"13gene" <- ifelse(newruleslhs$"13" %in% coop_vector,1,0) + ifelse(newruleslhs$"13" %in% self_vector,1,0)
  newruleslhs$"14gene" <- ifelse(newruleslhs$"14" %in% coop_vector,1,0) + ifelse(newruleslhs$"14" %in% self_vector,1,0)
  newruleslhs$"15gene" <- ifelse(newruleslhs$"15" %in% coop_vector,1,0) + ifelse(newruleslhs$"15" %in% self_vector,1,0)
  # newruleslhs$"16gene" <- ifelse(newruleslhs$"16" %in% coop_vector,1,0) + ifelse(newruleslhs$"16" %in% self_vector,1,0)
  # newruleslhs$"17gene" <- ifelse(newruleslhs$"17" %in% coop_vector,1,0) + ifelse(newruleslhs$"17" %in% self_vector,1,0)

  numcoop <- ifelse(newruleslhs$"1" %in% coop_vector,1,0) + ifelse(newruleslhs$"2" %in% coop_vector,1,0)+ ifelse(newruleslhs$"3" %in% coop_vector,1,0) + ifelse(newruleslhs$"4" %in% coop_vector,1,0)+ ifelse(newruleslhs$"5" %in% coop_vector,1,0)+ ifelse(newruleslhs$"6" %in% coop_vector,1,0)+ ifelse(newruleslhs$"7" %in% coop_vector,1,0)+ ifelse(newruleslhs$"8" %in% coop_vector,1,0)+ ifelse(newruleslhs$"9" %in% coop_vector,1,0)+ ifelse(newruleslhs$"10" %in% coop_vector,1,0)+ ifelse(newruleslhs$"11" %in% coop_vector,1,0)+ ifelse(newruleslhs$"12" %in% coop_vector,1,0)+ ifelse(newruleslhs$"13" %in% coop_vector,1,0)+ ifelse(newruleslhs$"14" %in% coop_vector,1,0)+ ifelse(newruleslhs$"15" %in% coop_vector,1,0)
  newruleslhs$numcoop <- numcoop
  rules$numcoop <- numcoop
  
  numself <- ifelse(newruleslhs$"1" %in% self_vector,1,0) + ifelse(newruleslhs$"2" %in% self_vector,1,0)+ ifelse(newruleslhs$"3" %in% self_vector,1,0) + ifelse(newruleslhs$"4" %in% self_vector,1,0)+ ifelse(newruleslhs$"5" %in% self_vector,1,0)+ ifelse(newruleslhs$"6" %in% self_vector,1,0)+ ifelse(newruleslhs$"7" %in% self_vector,1,0)+ ifelse(newruleslhs$"8" %in% self_vector,1,0)+ ifelse(newruleslhs$"9" %in% self_vector,1,0)+ ifelse(newruleslhs$"10" %in% self_vector,1,0)+ ifelse(newruleslhs$"11" %in% self_vector,1,0)+ ifelse(newruleslhs$"12" %in% self_vector,1,0)+ ifelse(newruleslhs$"13" %in% self_vector,1,0)+ ifelse(newruleslhs$"14" %in% self_vector,1,0)+ ifelse(newruleslhs$"15" %in% self_vector,1,0)
  newruleslhs$numself <- numself
  rules$numself <- numself
}

rules$pctself <- (rules$numself+1)/rules$size

saveRDS(rules,"../mcr_lr_model_oddsratios.RDS")
write.csv(rules %>% select(-genevector), "../mcr_lr_model_oddsratios.csv", row.names = FALSE)


###############################
## Coop/Self Ratio for MCR vs. Not
## Distribution of coopertivity of MCR-containing isolates
## Look in LR script to line 51

## Cols w/ coop-prefix/everything else (for mcr containing vs. not)
library(jsonlite)
library(stringr)
#library(glmnet)
library(caret)
#library(arules)

# Run these *nix commands
#sed 's/'\''//g' PDG000000004.1024.reference_target.tree.newick | sed 's/:-*[0-9]\.*[0-9]*\(e-[0-9]*\)*//g' | sed 's/,'\('/'\('/g' | sed 's/,/ /g' > e.coli.paren
#sed 's/'\)'*'\('*'\ '*PDT/,PDT/g' e.coli.paren | sed 's/^,//g' | sed 's/'\)'*;$//g' > e.coli.ids
#---------
pdt <- unique(names(read.csv("../e.coli.ids")))

#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=1000000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Casc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=asc)")
#saveRDS(json.stream, file = "e.coli.RDS")
json.stream <- readRDS("../e.coli.RDS")

id <- substr(json.stream[["ngout"]][["data"]][["content"]][["id"]], 19, 33)
creation_date_time <- as.POSIXct(json.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")
collection_year <- as.numeric(substr(json.stream[["ngout"]][["data"]][["content"]][["collection_date"]], 1, 4))
collection_year[is.na(collection_year)] <- as.numeric(format(creation_date_time[is.na(collection_year)], "%Y"))
location <- json.stream[["ngout"]][["data"]][["content"]][["geo_loc_name"]]
isolation_type <- json.stream[["ngout"]][["data"]][["content"]][["epi_type"]]
isolation_source <- json.stream[["ngout"]][["data"]][["content"]][["isolation_source"]]
fdf <- data.frame(id, creation_date_time, collection_year, location, isolation_type, isolation_source)
rdf <- fdf[fdf$id %in% pdt,]
fgenotypes <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(fgenotypes) <- id
rgenotypes <- fgenotypes
rgenotypes[!(names(fgenotypes) %in% pdt)] <- NULL
#fgenotypes[fgenotypes == "NULL"] <- NULL
#rgenotypes[rgenotypes == "NULL"] <- NULL
genotypes <- fgenotypes[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]
genoid <- id[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]

# ---------
# Genotypes
# ---------
genotypesall <- unique(sort(unlist(genotypes)))
genobinv <- lapply(genotypes, function(x){(genotypesall %in% x) * 1})
gdf <- as.data.frame(matrix(unlist(genobinv), ncol = length(genotypesall), byrow = T))
names(gdf) <- genotypesall
gdf.mcr <- gdf[,-grep("mcr", names(gdf))]
gdf.mcr$mcr <- (rowSums(gdf[,grep("mcr", names(gdf))]) > 0) * 1
# gdf.mcr$mcr <- as.factor(gdf.mcr$mcr)

allisolates <- gdf.mcr
nonmcrisolates <- gdf.mcr %>% filter(mcr == 0)
onlymcrisolates <- gdf.mcr %>% filter(mcr == 1)
# onlyblaisolates <- gdf.mcr %>% filter_at(vars(starts_with("bla")), any_vars(. == 1))
# onlyrmtisolates <- gdf.mcr %>% filter_at(vars(starts_with("rmt")), any_vars(. == 1))
# onlyermisolates <- gdf.mcr %>% filter_at(vars(starts_with("erm")), any_vars(. == 1))


## Identify genes are their respective socialities
gene <- read.csv("../AMR_FunctionalMechanisms.csv", header=T)

for (i in 1:nrow(gene)) {
  coop = ifelse(gene$Sociality == "Cooperative", paste(gene$AMR.Gene), NA)
  self = ifelse(gene$Sociality == "Selfish", paste(gene$AMR.Gene), NA)
  unkn = ifelse(gene$Sociality == "", paste(gene$AMR.Gene), NA)
  gene$coop <- coop
  gene$self <- self
  gene$unkn <- unkn
}

## Makes vectors contaning genes in their categories
coop_vector <- c(gene$coop[!is.na(gene$coop)])
self_vector <- c(gene$self[!is.na(gene$self)],"mcr")
unkn_vector <- c(gene$self[!is.na(gene$unkn)])

coop_vector_ors <- paste(coop_vector, collapse = "|")
self_vector_ors <- paste(self_vector, collapse = "|")
unkn_vector_ors <- paste(unkn_vector, collapse = "|")
all_vector_ors <- paste0(coop_vector_ors,"|",self_vector_ors,"|",unkn_vector_ors)

## Get isolates with at least 1 individualistic gene
onlyselfisolates <- gdf.mcr %>% filter_at(vars(matches(self_vector_ors)), any_vars(. == 1))

## Look at ratios overall vs. only MCR vs only ...

allisolates <- allisolates %>% 
  mutate(set = "all",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

hist(allisolates$pctself)


onlyselfisolates <- onlyselfisolates %>% 
  mutate(set = "self",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

hist(onlyselfisolates$pctself)

nonmcrisolates <- nonmcrisolates %>% 
  mutate(set = "non-mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

hist(nonmcrisolates$pctself)

onlymcrisolates <- onlymcrisolates %>% 
  mutate(set = "mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

hist(onlymcrisolates$pctself)


# onlyblaisolates <- onlyblaisolates %>% 
#   mutate(set = "bla",
#          count = rowSums(select(., matches(all_vector_ors))),
#          numcoop = rowSums(select(., matches(coop_vector_ors))),
#          numself = rowSums(select(., matches(self_vector_ors))),
#          numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
#   mutate(ratio = numcoop/(numself),
#          pctcoop = numcoop/count,
#          pctself = numself/count)
# 
# hist(onlyblaisolates$pctself)
# 
# 
# onlyrmtisolates <- onlyrmtisolates %>% 
#   mutate(set = "rmt",
#          count = rowSums(select(., matches(all_vector_ors))),
#          numcoop = rowSums(select(., matches(coop_vector_ors))),
#          numself = rowSums(select(., matches(self_vector_ors))),
#          numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
#   mutate(ratio = numcoop/(numself),
#          pctcoop = numcoop/count,
#          pctself = numself/count)
# 
# hist(onlyrmtisolates$pctself)
# 
# onlyermisolates <- onlyermisolates %>% 
#   mutate(set = "erm",
#          count = rowSums(select(., matches(all_vector_ors))),
#          numcoop = rowSums(select(., matches(coop_vector_ors))),
#          numself = rowSums(select(., matches(self_vector_ors))),
#          numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
#   mutate(ratio = numcoop/(numself),
#          pctcoop = numcoop/count,
#          pctself = numself/count)
# 
# hist(onlyermisolates$pctself)


## Analyze Distributions of ALL vs. MCR
# devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)
library(ggsci)
library(ggplot2)

mcr_vs_nonmcr_self_vs_all_plotdata <- rbind(allisolates, onlyselfisolates, nonmcrisolates, onlymcrisolates)

write.csv(mcr_vs_nonmcr_self_vs_all_plotdata,
          file = "mcr_vs_nonmcr_self_vs_all_plotdata.csv")

mcr_vs_nonmcr_self_vs_all_plotdata <- read.csv("mcr_vs_nonmcr_self_vs_all_plotdata.csv")

ggplot2.histogram(data=mcr_vs_nonmcr_self_vs_all_plotdata,
                  xName='pctself',
                  groupName='set',
                  xtitle = "Percentage of Selfish Genes",
                  ytitle = "Density",
                  legendPosition="top",
                  alpha=0.75,
                  addDensity=TRUE) + 
  scale_color_rickandmorty() + 
  scale_fill_rickandmorty()

### Stats tests
shapiro.test(sample(mcr_vs_nonmcr_self_vs_all_plotdata$pctself,size = 5000))
## (NOT NORMAL)
mcr_vs_nonmcr_self_vs_all_plotdata$set <- factor(mcr_vs_nonmcr_self_vs_all_plotdata$set) 

### KW Tests
## All
kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata)
## MCR vs Non-MCR
kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("self", "all")))
## Selfish vs. All
kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("mcr", "non-mcr")))
## MCR vs Selfish
kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("non-mcr", "all")))
## Non-MCR vs Selfish
kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("mcr", "all")))


### MWU Tests

## Selfish vs. All
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("mcr", "non-mcr")),
            paired = FALSE)

## MCR vs. All
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("self", "non-mcr")),
            paired = FALSE)

## Non-MCR vs. All
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("self", "mcr")),
            paired = FALSE)

## MCR vs Non-MCR
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("self", "all")),
            paired = FALSE)

## MCR vs Selfish
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("non-mcr", "all")),
            paired = FALSE)

## Non-MCR vs Selfish
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(!set %in% c("mcr", "all")),
            paired = FALSE)
