#######################
## Script for finding all PDT numbers for publication crawling
## By: Colby T. Ford, Ph.D., Gabriel Zenarosa, Ph.D., David Brown, Kevin Smith, and Daniel Janies, Ph.D.
#######################

## Load in Packages
library(jsonlite)
library(stringr)
library(dplyr)
library(tidyr)
library(arules)

## Generate New Extract from the NCBI Database
#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=100000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Cdesc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=desc)")
#saveRDS(json.stream, file = "json.stream.RDS") ## to 1/8/2019

## Read in Data
json.stream <- readRDS("../e.coli.RDS")

amr.list.raw <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(amr.list.raw) <- json.stream[["ngout"]][["data"]][["content"]][["id"]]
creation_date_time <- as.POSIXct(json.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")

## Filter to cases after 4/4/2016
amr.list.raw <- amr.list.raw[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]

amr.list.raw[amr.list.raw == "NULL"] <- NULL
amr.all <- sort(unique(unlist(amr.list.raw)))

## Find all mcr variants
# mcr.variants <- amr.all[str_detect(amr.all, "mcr-1|mcr-2")] %>% 
mcr.variants <- amr.all[str_detect(amr.all, "mcr-")] %>% 
  sort(decreasing = TRUE) %>% 
  paste0(collapse = "|")

## Convert all mcr genotype variants to "mcr"
mcr.subset <- amr.list.raw %>% 
  lapply(., str_replace_all, pattern = mcr.variants, replacement = "mcr") %>% 
  lapply(., unique) #Removes duplicate "mcr" genotypes from previous step



## Read in Rule Sets
ARM_rules <- readRDS("../mcr_rules_all.RDS") %>% 
  mutate(method = "arm",
         id = row.names(.)) %>% 
  select(-rule)

LR_rules <- readRDS("../mcr_lr_model_oddsratios.RDS") %>%
  mutate(lhs = stringr::str_replace_all(GeneSet,":",","),
         rhs = "mcr",
         id = row.names(.))

ADASYN_rules <- LR_rules %>% 
  mutate(oddsratio = OddsRatio_adasyn,
         method = "adasyn") %>% 
  select(-OddsRatio_glm, -OddsRatio_adasyn, -genevector, -GeneSet, -numcoop, -numself, -pctself) %>% 
  na.omit()

GLM_rules <- LR_rules %>% 
  mutate(oddsratio = OddsRatio_glm,
         method = "glm") %>% 
  select(-OddsRatio_glm, -OddsRatio_adasyn, -genevector, -GeneSet, -numcoop, -numself, -pctself) %>% 
  na.omit()


mcr_rules_df <- bind_rows(ARM_rules, ADASYN_rules, GLM_rules)

readr::write_csv(mcr_rules_df, "all_rules.csv")

## Create strings of each genotype
mcrgenotypes <- sapply(mcr.subset, function(x){unlist(x) %>% paste0(collapse = ",")}) %>%
  as.data.frame()
colnames(mcrgenotypes) <- "genotype"

## Create Empty Dataframe of All PDTs with MCR
matches <- data.frame(pdt = names(mcr.subset),
                      genotype = mcrgenotypes$genotype
)

matches <- data.frame(method = character(0),
                      id = numeric(0),
                      pdt = character(0),
                      genotype = character(0))

## Create a function to compare the rule vs. the genotype sets
comparefxn <- function(a, b){
  if(setequal(a, unlist(b))){
    response <- "Match"
  } else if(all(a %in% b)){
    response <- "Subset"
  } else if(all(b %in% a)){
    response <- "Superset"
  } else{
    response <- NA
  }
  return(response)
}

## Loop through each rule and see if any exist in the MCR Subset
for (i in 1:nrow(mcr_rules_df)){
  rulename <- paste0(mcr_rules_df$lhs[i], "=>", mcr_rules_df$rhs[i])
  rulelhs <- mcr_rules_df$lhs[i] %>% 
    strsplit(split=",") %>% 
    unlist()
  
  cat(i, ". ", rulename, "\n", sep = "")
  
  # itermatches <- sapply(mcr.subset, function(x){setequal(rulelhs, unlist(x))}) %>% 
  #   as.data.frame()
  
  itermatches <- sapply(mcr.subset, function(x){comparefxn(rulelhs, x)}) %>% 
    as.data.frame()
  
  colnames(itermatches) <- rulename
  rownames(itermatches) <- names(mcr.subset)
  matches <- cbind(matches, itermatches) ## Append to dataframe
}

# saveRDS(matches, "matches.RDS")
# matches <- readRDS("matches.RDS")

## Check to see if there are any matches
# for (i in 2:ncol(matches)){
#   cat(i,
#       ". ",
#       colnames(matches[,i]),
#       ": ",
#       print(levels(as.factor(matches[,i]))),
#       "\n",
#       sep = "")
# } 

## Reshape to Get only Subsets/Matches
library(reshape2)

outputpdts <- melt(matches, id = c("pdt", "genotype")) %>% 
  filter(value %in% c("Match","Subset")) %>% 
  select(-one_of("value"))

colnames(outputpdts) <- c("id", "genotype", "mcr_group")

## Separate IDs into PDG and PDT numbers
outputpdts <- outputpdts %>%
  tidyr::separate(id, c("pdg", "pdt"), sep = "_") %>%
  select(-c("pdg","genotype"))

readr::write_csv(outputpdts, "allPDTs.csv")
