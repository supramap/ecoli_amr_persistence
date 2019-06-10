############################################################
## Script for generating Association Rules for            ##
## E. coli genotype sets containing `mcr`                 ##
## By: Colby T. Ford, Ph.D., Gabriel Zenarosa, Ph.D.,     ##
##     David Brown, Kevin Smith, and Daniel Janies, Ph.D. ##
############################################################

## Load in Packages
library(jsonlite)
library(stringr)
library(dplyr)
library(tidyr)
library(arules)

## Generate New Extract from the NCBI Database
#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=100000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Cdesc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=desc)")
#saveRDS(json.stream, file = "json.stream.RDS")

## Read in Data
json.stream <- readRDS("e.coli.RDS")

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

## Filter to strains containing "mcr"
#mcr.subset <- mcr.subset[unlist(lapply(mcr.subset, function(x){"mcr" %in% x}))]

genotypes <- sort(unique(unlist(mcr.subset)))

## Convert to transactions for the a priori algorithm
mcr.transactions <- as(mcr.subset, 'transactions')


######################################
## Generate Rules with the pattern {*|-mcr} => {mcr}
# 
# mcr_rules<- apriori(mcr.transactions,
#                     parameter = list(minlen = 2,
#                                      maxlen = 14,
#                                      support = (200/length(genotypes))-1,
#                                      confidence =  0.5),
#                     appearance = list(default = "none",
#                                       lhs = genotypes[which(genotypes != "mcr")],
#                                       rhs = "mcr"),
#                     control = list(memopt = FALSE)
# )
# 
# mcr_rules_df <- data.frame(
#   lhs = labels(lhs(mcr_rules)),
#   rhs = labels(rhs(mcr_rules)),
#   mcr_rules@quality
# )
# 
# # Remove brackets and convert all sets to lists
# mcr_rules_df$lhs <- mcr_rules_df$lhs %>% 
#   str_replace_all(c("\\{|\\}"),"")# %>%
#   #strsplit(split = ",")
# 
# mcr_rules_df$rhs <- mcr_rules_df$rhs %>% 
#   str_replace_all(c("\\{|\\}"),"")
# 
# ## Flag which rules are supersets
# # mcr_rules_df$is_superset <- FALSE
# # 
# # for (i in 1:length(mcr_rules_df$lhs)){
# #   for (j in 1:length(mcr_rules_df$lhs)){
# #     cat("Checking row",i,"against row",j,"\n")
# #     #if (all(mcr_rules_df$lhs[j] %in% mcr_rules_df$lhs[i])){
# #     if (setequal(intersect(mcr_rules_df$lhs[j],
# #                            mcr_rules_df$lhs[i]),
# #                  mcr_rules_df$lhs[i])){
# #       mcr_rules_df$is_superset[i] <- TRUE
# #     }
# #   }
# # }
# 
# ## Write out the rules results
# saveRDS(mcr_rules_df,
#         file = "mcr_rules.RDS")
# 
# write.csv(mcr_rules_df,
#           file = "mcr_rules.csv",
#           append = FALSE)

###########################################
## Generate Rules Iteratively (with the pattern {*|-mcr} => {mcr})

## Define rule sizes and isolate counts to explore
rulesizes <- seq(2, 14, by = 1)
isolatesizes <- lapply(amr.list.raw, length) %>% 
  matrix() %>% 
  data.frame()

colnames(isolatesizes) <- c("size")
isolatesizes$size <- as.factor(as.character(isolatesizes$size))

isolatesizes <- isolatesizes %>%
  group_by(size) %>% 
  tally()

isolatesizes$size <- as.numeric(isolatesizes$size)

isolatesizes <- isolatesizes %>%
  arrange() %>% 
  filter(size <= max(rulesizes) & size >= min(rulesizes))

rulespace <- data.frame(rulesize = rulesizes,
                        isolatesize = isolatesizes$n)
# rulesizes <- c(2,4,5,8,10,14)
# isolatesizes <- c(1944, 436, 265, 183, 11, 11)
# rulespace <- data.frame(rulesize = rulesizes,
#                         isolatesize = isolatesizes)

## Create empty dataframe
mcr_rules_df <- data.frame(
  lhs = character(0),
  rhs = character(0),
  size = numeric(0),
  support = numeric(0),
  confidence = numeric(0),
  lift = numeric(0),
  count = numeric(0)
)

for (i in 1:nrow(rulespace)){
  ## Get rule size and support metric for this iteration
  size <- rulespace$rulesize[i]
  numerator <- rulespace$isolatesize[i]
  
  cat("Testing rules of size:", size, "\n")
  
  ## Generate MCR Rules
  mcr_rules<- apriori(mcr.transactions,
                      parameter = list(minlen = size,
                                       maxlen = size,
                                       #support = 0.002,
                                       #support = 5*(numerator/length(amr.list.raw)),
                                       support = 0.5*(numerator/length(amr.list.raw)),
                                       confidence =  2/3,
                                       maxtime = 60),
                      appearance = list(default = "none",
                                        lhs = genotypes[which(genotypes != "mcr")],
                                        rhs = "mcr"),
                      control = list(memopt = FALSE)
  )
  
  if (length(mcr_rules) > 0){
    ## Convert to dataframe
    mcr_rules_df_iter <- data.frame(
      lhs = labels(lhs(mcr_rules)),
      rhs = labels(rhs(mcr_rules)),
      size = size,
      mcr_rules@quality
    )
    
    ## Append this iterations's results to main dataframe
    mcr_rules_df <- rbind(mcr_rules_df,
                          mcr_rules_df_iter)
  }
}



## Remove brackets and convert all sets to lists
mcr_rules_df$lhs <- mcr_rules_df$lhs %>% 
  str_replace_all(c("\\{|\\}"),"") #%>%
  #strsplit(split=",")

mcr_rules_df$rhs <- mcr_rules_df$rhs %>% 
  str_replace_all(c("\\{|\\}"),"")

mcr_rules_df$rule <- paste0(mcr_rules_df$lhs, ",", mcr_rules_df$rhs) %>%
  strsplit(split=",")

## Write out the rules results
saveRDS(mcr_rules_df,
        file = "mcr_rules.RDS")

write.csv(mcr_rules_df,
          file = "mcr_rules.csv",
          append = FALSE)


#############################
## Get PDTs that match rules

## Create strings of each genotype
mcrgenotypes <- sapply(mcr.subset, function(x){unlist(x) %>% paste0(collapse = ",")}) %>%
  as.data.frame()
colnames(mcrgenotypes) <- "genotype"

## Create Empty Dataframe of All PDTs with MCR
matches <- data.frame(pdt = names(mcr.subset),
                      genotype = mcrgenotypes$genotype)

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

## Check to see if there are any matches
for (i in 2:ncol(matches)){
  cat(i,
      ". ",
      colnames(matches[,i]),
      ": ",
      print(levels(as.factor(matches[,i]))),
      "\n",
      sep = "")
} 

## Reshape to Get only Subsets/Supersets/Matches
library(reshape2)

outputpdts <- melt(matches, id = c("pdt", "genotype")) %>% 
  filter(value == "Subset") %>% 
  select(-one_of("value"))

colnames(outputpdts) <- c("id", "genotype", "mcr_group")

## Separate IDs into PDG and PDT numbers
outputpdts <- outputpdts %>%
  separate(id, c("pdg", "pdt"), sep = "_") %>%
  select(-c("pdg","genotype"))

readr::write_csv(outputpdts, "MCR_Groups_ARM_PDTs.csv")

###########################################
## Validation Set - Create rules using new data to compare

## New data From Jan 9, 2019 to May 8, 2019
json.stream <- readRDS("e.coli_validation.RDS")

amr.list.raw <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(amr.list.raw) <- json.stream[["ngout"]][["data"]][["content"]][["id"]]
creation_date_time <- as.POSIXct(json.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")

## Filter to cases after 1/9/2019
amr.list.raw <- amr.list.raw[creation_date_time >= as.POSIXct("2019-01-09 00:00:00", format = "%Y-%m-%d %H:%M:%S")]

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

## Filter to strains containing "mcr"
#mcr.subset <- mcr.subset[unlist(lapply(mcr.subset, function(x){"mcr" %in% x}))]

genotypes <- sort(unique(unlist(mcr.subset)))

## Convert to transactions for the a priori algorithm
mcr.transactions <- as(mcr.subset, 'transactions')

## Generate Rules Iteratively (with the pattern {*|-mcr} => {mcr})
## Define rule sizes and isolate counts to explore
rulesizes <- seq(2, 14, by = 1)
isolatesizes <- lapply(amr.list.raw, length) %>% 
  matrix() %>% 
  data.frame()

colnames(isolatesizes) <- c("size")
isolatesizes$size <- as.factor(as.character(isolatesizes$size))

isolatesizes <- isolatesizes %>%
  group_by(size) %>% 
  tally()

isolatesizes$size <- as.numeric(isolatesizes$size)

isolatesizes <- isolatesizes %>%
  arrange() %>% 
  filter(size <= max(rulesizes) & size >= min(rulesizes))

rulespace <- data.frame(rulesize = rulesizes,
                        isolatesize = isolatesizes$n)
# rulesizes <- c(2,4,5,8,10,14)
# isolatesizes <- c(1944, 436, 265, 183, 11, 11)
# rulespace <- data.frame(rulesize = rulesizes,
#                         isolatesize = isolatesizes)

## Create empty dataframe
mcr_rules_df <- data.frame(
  lhs = character(0),
  rhs = character(0),
  size = numeric(0),
  support = numeric(0),
  confidence = numeric(0),
  lift = numeric(0),
  count = numeric(0)
)

for (i in 1:nrow(rulespace)){
  ## Get rule size and support metric for this iteration
  size <- rulespace$rulesize[i]
  numerator <- rulespace$isolatesize[i]
  
  cat("Testing rules of size:", size, "\n")
  
  ## Generate MCR Rules
  mcr_rules<- apriori(mcr.transactions,
                      parameter = list(minlen = size,
                                       maxlen = size,
                                       #support = 0.002,
                                       #support = 5*(numerator/length(amr.list.raw)),
                                       support = 0.5*(numerator/length(amr.list.raw)),
                                       confidence =  2/3,
                                       maxtime = 60),
                      appearance = list(default = "none",
                                        lhs = genotypes[which(genotypes != "mcr")],
                                        rhs = "mcr"),
                      control = list(memopt = FALSE)
  )
  
  if (length(mcr_rules) > 0){
    ## Convert to dataframe
    mcr_rules_df_iter <- data.frame(
      lhs = labels(lhs(mcr_rules)),
      rhs = labels(rhs(mcr_rules)),
      size = size,
      mcr_rules@quality
    )
    
    ## Append this iterations's results to main dataframe
    mcr_rules_df <- rbind(mcr_rules_df,
                          mcr_rules_df_iter)
  }
}


## Remove brackets and convert all sets to lists
mcr_rules_df$lhs <- mcr_rules_df$lhs %>% 
  str_replace_all(c("\\{|\\}"),"") #%>%
#strsplit(split=",")

mcr_rules_df$rhs <- mcr_rules_df$rhs %>% 
  str_replace_all(c("\\{|\\}"),"")

mcr_rules_df$rule <- paste0(mcr_rules_df$lhs, ",", mcr_rules_df$rhs) %>%
  strsplit(split=",")

## Write out the rules results
saveRDS(mcr_rules_df,
        file = "mcr_validation_rules.RDS")

write.csv(mcr_rules_df %>% select(-rule),
          file = "mcr_validation_rules.csv",
          append = FALSE)
