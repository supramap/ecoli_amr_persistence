#######################
## Script for generating Association Rules for E. coli genotype sets containing MCR
## Colby T. Ford, Ph.D.
#######################

## Load in Packages
library(jsonlite)
library(stringr)
library(dplyr)
library(tidyr)
library(arules)

## Generate New Extract from the NcBI Database
#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=100000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Cdesc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=desc)")
#saveRDS(json.stream, file = "json.stream.RDS")

## Read in Data
json.stream <- readRDS("e.coli.RDS")

amr.list.raw <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(amr.list.raw) <- json.stream[["ngout"]][["data"]][["content"]][["id"]]
amr.list.raw[amr.list.raw == "NULL"] <- NULL
amr.all <- sort(unique(unlist(amr.list.raw)))

## Find all MCR variants
mcr.variants <- amr.all[str_detect(amr.all, "mcr-1|mcr-2")] %>% 
  sort(decreasing = TRUE) %>% 
  paste0(collapse = "|")

## Convert all mcr genotype variants to "mcr"
mcr.subset <- amr.list.raw %>% 
  lapply(., str_replace_all, pattern = mcr.variants, replacement = "mcr") %>% 
  lapply(., unique)

## Filter strains containing "mcr"
mcr.subset <- mcr.subset[unlist(lapply(mcr.subset, function(x){"mcr" %in% x}))]

genotypes <- sort(unique(unlist(mcr.subset)))

## Convert to transactions for the a priori algorithm
mcr.transactions <- as(mcr.subset, 'transactions')

######################################
## Generate Rules with the pattern {*|-mcr} => {mcr}

mcr_rules<- apriori(mcr.transactions,
                    parameter = list(minlen = 2,
                                     maxlen = 14,
                                     support = 200.0/length(genotypes),
                                     confidence =  0.5),
                    appearance = list(default = "none",
                                      lhs = genotypes[which(genotypes != "mcr")],
                                      rhs = "mcr"),
                    control = list(memopt = FALSE)
)

mcr_rules_df <- data.frame(
  lhs = labels(lhs(mcr_rules)),
  rhs = labels(rhs(mcr_rules)),
  mcr_rules@quality
)

# Remove brackets and convert all sets to lists
mcr_rules_df$lhs <- mcr_rules_df$lhs %>% 
  str_replace_all(c("\\{|\\}"),"")# %>%
  #strsplit(split = ",")

mcr_rules_df$rhs <- mcr_rules_df$rhs %>% 
  str_replace_all(c("\\{|\\}"),"")

## Flag which rules are supersets
# mcr_rules_df$is_superset <- FALSE
# 
# for (i in 1:length(mcr_rules_df$lhs)){
#   for (j in 1:length(mcr_rules_df$lhs)){
#     cat("Checking row",i,"against row",j,"\n")
#     #if (all(mcr_rules_df$lhs[j] %in% mcr_rules_df$lhs[i])){
#     if (setequal(intersect(mcr_rules_df$lhs[j],
#                            mcr_rules_df$lhs[i]),
#                  mcr_rules_df$lhs[i])){
#       mcr_rules_df$is_superset[i] <- TRUE
#     }
#   }
# }

## Write out the rules results
saveRDS(mcr_rules_df,
        file = "mcr_rules.RDS")

write.csv(mcr_rules_df,
          file = "mcr_rules.csv",
          append = FALSE)


