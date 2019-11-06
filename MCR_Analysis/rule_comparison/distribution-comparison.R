############################################################
## Scripts for comparing sociation distribution           ##
## between different isolate subsets                      ##
## By:  Kevin Smith, Colby T. Ford, Ph.D.,                ##
##      Gabriel Zenarosa, Ph.D., David Brown,             ##
##      and Daniel Janies, Ph.D.                          ##
############################################################

###############################
## Coop/Self Ratio for MCR vs. Not
## Distribution of coopertivity of MCR-containing isolates

## Cols w/ coop-prefix/everything else (for mcr containing vs. not)
library(jsonlite)
library(stringr)
library(dplyr)
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

collection_year_filtered <- collection_year[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]
location_filtered <- location[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]
isolation_type_filtered <- isolation_type[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")]
isolation_source_filtered <- isolation_source[creation_date_time >= as.POSIXct("2016-04-04 20:14:38", format = "%Y-%m-%d %H:%M:%S")] 

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

## Subset Data for Distribution Analyses
# allisolates <- gdf.mcr
allisolates <- data.frame(id = genoid,
                          collection_year = collection_year_filtered,
                          location = location_filtered,
                          isolation_type = isolation_type_filtered,
                          isolation_source = isolation_source_filtered,
                          gdf.mcr)
allisolates <- allisolates %>% 
  tidyr::separate(location,
                  c("country", NA),
                  sep = ":",
                  extra = "drop")

## By Country
chinaisolates <- allisolates %>% filter(country == "China")
nonchinaisolates <- allisolates %>% filter(country != "China")

## By Isolate Type
clinicalisolates <- allisolates %>% filter(isolation_type == "clinical")
environmentalisolates <- allisolates %>% filter(isolation_type == "environmental/other")
clinicalmcrisolates <- allisolates %>% filter(isolation_type == "clinical") %>% filter(mcr == 1)
environmentalmcrisolates <- allisolates %>% filter(isolation_type == "environmental/other") %>% filter(mcr == 1)

## By Gene Type Subset
nonmcrisolates <- allisolates %>% filter(mcr == 0)
onlymcrisolates <- allisolates %>% filter(mcr == 1)

## Identify genes are their respective socialities
gene <- read.csv("../AMR_FunctionalMechanisms.csv", header=T) %>% 
  filter(Gene.Category != "mcr") # To avoid counting mcr as individualistic

mcrs <- read.csv("../AMR_FunctionalMechanisms.csv", header=T) %>% 
  filter(Gene.Category == "mcr")

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
# self_vector <- c(gene$self[!is.na(gene$self)],"mcr") # To count mcr as individualistic
self_vector <- c(gene$self[!is.na(gene$self)])
unkn_vector <- c(gene$unkn[!is.na(gene$unkn)])
mcr_vector <- c(as.character(mcrs$AMR.Gene))

coop_vector_ors <- paste(coop_vector, collapse = "|")
self_vector_ors <- paste(self_vector, collapse = "|")
unkn_vector_ors <- paste(unique(unkn_vector), collapse = "|")
mcr_vector_ors <- paste(mcr_vector, collapse = "|")
all_vector_ors <- paste0(coop_vector_ors,"|",self_vector_ors,"|",unkn_vector_ors,"|",mcr_vector_ors)

## Get isolates with at least 1 individualistic gene
onlyselfisolates <- allisolates %>% filter_at(vars(matches(self_vector_ors)), any_vars(. == 1))
selfwithoutmcr <- onlyselfisolates %>% filter(mcr == 0)

## Exclude counting MCR as Individualistic
mcrisolates_notcountingmcr <- allisolates %>% filter(mcr == 1) %>% mutate(mcr = 0)
onlyselfisolates_notcountingmcr <- allisolates %>% mutate(mcr = 0) %>% filter_at(vars(matches(self_vector_ors)), any_vars(. == 1))

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
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(allisolates$pctself)

## By Country
chinaisolates <- chinaisolates %>% 
  mutate(set = "chineseisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(chinaisolates$pctself)

nonchinaisolates <- nonchinaisolates %>% 
  mutate(set = "nonchineseisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(nonchinaisolates$pctself)

## By Isolate Type
clinicalisolates <- clinicalisolates %>% 
  mutate(set = "clinicalisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(clinicalisolates$pctself)

environmentalisolates <- environmentalisolates %>% 
  mutate(set = "environmentalisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(environmentalisolates$pctself)

clinicalmcrisolates <- clinicalmcrisolates %>% 
  mutate(set = "clinicalmcrisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(clinicalmcrisolates$pctself)

environmentalmcrisolates <- environmentalmcrisolates %>% 
  mutate(set = "environmentalmcrisolates",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(environmentalmcrisolates$pctself)

## MCR vs Non-MCR vs Individualistic vs All
onlyselfisolates <- onlyselfisolates %>% 
  mutate(set = "self",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(onlyselfisolates$pctself)

nonmcrisolates <- nonmcrisolates %>% 
  mutate(set = "non-mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(nonmcrisolates$pctself)

onlymcrisolates <- onlymcrisolates %>% 
  mutate(set = "mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(onlymcrisolates$pctself)

selfwithoutmcr <- selfwithoutmcr %>% 
  mutate(set = "self_without_mcr",
         count = rowSums(select(., matches(all_vector_ors))),
         numcoop = rowSums(select(., matches(coop_vector_ors))),
         numself = rowSums(select(., matches(self_vector_ors))),
         numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  mutate(ratio = numcoop/(numself),
         pctcoop = numcoop/count,
         pctself = numself/count) %>% 
  select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)

# hist(selfwithoutmcr$pctself)

## ## MCR vs Non-MCR vs Individualistic vs All (not counting MCR as Individualistic)
# mcrisolates_notcountingmcr <- mcrisolates_notcountingmcr %>% 
#   mutate(set = "mcr_notcountingmcr",
#          count = rowSums(select(., matches(all_vector_ors))),
#          numcoop = rowSums(select(., matches(coop_vector_ors))),
#          numself = rowSums(select(., matches(self_vector_ors))),
#          numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
#   mutate(ratio = numcoop/(numself),
#          pctcoop = numcoop/count,
#          pctself = numself/count) %>% 
#   select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)
# 
# hist(mcrisolates_notcountingmcr$pctself)
# 
# onlyselfisolates_notcountingmcr <- onlyselfisolates_notcountingmcr %>% 
#   mutate(set = "self_without_mcr_notcountingmcr",
#          count = rowSums(select(., matches(all_vector_ors))),
#          numcoop = rowSums(select(., matches(coop_vector_ors))),
#          numself = rowSums(select(., matches(self_vector_ors))),
#          numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
#   mutate(ratio = numcoop/(numself),
#          pctcoop = numcoop/count,
#          pctself = numself/count) %>% 
#   select(id, collection_year, country, isolation_type, isolation_source, set, count, numcoop, numself, numunkn, ratio, pctcoop, pctself)
# 
# hist(onlyselfisolates_notcountingmcr$pctself)

#######################
## Analyze Distributions of ALL vs. MCR

## MCR vs Non-MCR vs Individualistic vs All
mcr_vs_nonmcr_self_vs_all_plotdata <- rbind(allisolates,
                                            onlyselfisolates,
                                            nonmcrisolates,
                                            onlymcrisolates,
                                            selfwithoutmcr) # %>% 
# mutate(country = stringr::str_extract(location, "^[A-za-z ]+"))

write.csv(mcr_vs_nonmcr_self_vs_all_plotdata,
          file = "mcr_vs_nonmcr_self_vs_all_plotdata.csv")

mcr_vs_nonmcr_self_vs_all_plotdata <- read.csv("mcr_vs_nonmcr_self_vs_all_plotdata.csv")

## By Country
chinese_vs_non_chinese_plotdata <- rbind(chinaisolates,
                                         nonchinaisolates)

write.csv(chinese_vs_non_chinese_plotdata,
          file = "chinese_vs_non_chinese_plotdata.csv")

## By Isolation Type
clinical_vs_environmental_plotdata <- rbind(clinicalisolates,
                                            environmentalisolates)

write.csv(clinical_vs_environmental_plotdata,
          file = "clinical_vs_environmental_plotdata.csv")

clinical_vs_environmental_MCR_plotdata <- rbind(clinicalmcrisolates,
                                                environmentalmcrisolates)

write.csv(clinical_vs_environmental_MCR_plotdata,
          file = "clinical_vs_environmental_MCR_plotdata.csv")

## MCR vs Non-MCR vs Individualistic vs All (not counting MCR as Individualistic)
# mcr_vs_nonmcr_vs_self_notcountingmcr_plotdata <- rbind(mcrisolates_notcountingmcr,
#                                                        selfwithoutmcr_notcountingmcr,
#                                                        nonmcrisolates)
# 
# write.csv(mcr_vs_nonmcr_vs_self_notcountingmcr_plotdata,
#           file = "mcr_vs_nonmcr_vs_self_notcountingmcr_plotdata.csv")

# ### Stats tests
# shapiro.test(sample(mcr_vs_nonmcr_self_vs_all_plotdata$pctself,size = 5000))
# ## (NOT NORMAL)
# mcr_vs_nonmcr_self_vs_all_plotdata$set <- factor(mcr_vs_nonmcr_self_vs_all_plotdata$set) 

# ### KW Tests
# ## All
# kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata)
# ## MCR vs Non-MCR
# kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("self", "all")))
# ## Selfish vs. All
# kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("mcr", "non-mcr")))
# ## MCR vs Selfish
# kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("non-mcr", "all")))
# ## Non-MCR vs Selfish
# kruskal.test(pctself ~ set, data = mcr_vs_nonmcr_self_vs_all_plotdata %>% filter(!set %in% c("mcr", "all")))


### MWU Tests
## MCR vs Non-MCR
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "non-mcr")),
            paired = FALSE)


## MCR vs Non-MCR_Selfish
wilcox.test(pctself ~ set,
            data = mcr_vs_nonmcr_self_vs_all_plotdata %>%
              filter(set %in% c("mcr", "self_without_mcr")),
            paired = FALSE)

## Chinese vs. Non-Chinese Isolates
wilcox.test(pctself ~ set,
            data = chinese_vs_non_chinese_plotdata %>%
              filter(set %in% c("chineseisolates", "nonchineseisolates")),
            paired = FALSE)

## Clinical vs. Environmental Isolates
wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_plotdata,
            paired = FALSE)

## Clinical vs. Environmental MCR Isolates
wilcox.test(pctself ~ set,
            data = clinical_vs_environmental_MCR_plotdata,
            paired = FALSE)

# ## MCR vs Non-MCR (not counting MCR as Individualistic)
# wilcox.test(pctself ~ set,
#             data = mcr_vs_nonmcr_vs_self_notcountingmcr_plotdata %>% 
#               filter(set %in% c("mcr_notcountingmcr", "non-mcr")),
#             paired = FALSE)
# 
# ## MCR vs Non-MCR_Selfish (not counting MCR as Individualistic)
# wilcox.test(pctself ~ set,
#             data = mcr_vs_nonmcr_vs_self_notcountingmcr_plotdata %>% 
#               filter(set %in% c("mcr_notcountingmcr", "self_without_mcr_notcountingmcr")),
#             paired = FALSE)

################
## Subset Comparisons for Isolates with Gene Sets Matching Rules

rules <- read.csv("../web_crawler/all_rules.csv") %>% 
  mutate(lhs = stringr::str_replace_all(lhs, "[^A-Za-z0-9,]", ".")) %>% 
  mutate(lhs = stringr::str_split(lhs, pattern=","))

allisolates_sets <- allisolates %>% 
  # mutate(set = "all",
  #        count = rowSums(select(., matches(all_vector_ors))),
  #        numcoop = rowSums(select(., matches(coop_vector_ors))),
  #        numself = rowSums(select(., matches(self_vector_ors))),
  #        numunkn = rowSums(select(., matches(unkn_vector_ors)))) %>% 
  # mutate(ratio = numcoop/(numself),
  #        pctcoop = numcoop/count,
  #        pctself = numself/count) %>% 
  select(-one_of(c("collection_year",
                   "country",
                   "isolation_type",
                   "isolation_source"))) %>% 
  mutate(set = NA)

allisolates_sets$set <- genotypes
allisolates_sets$mcr <- gdf.mcr$mcr

genes <- colnames(allisolates_sets %>% select(-one_of(c("id","set","count","numcoop","numself","numunkn","ratio","pctcoop","pctself"))))

# for (i in 1:nrow(allisolates_sets)){
#   cat(paste0(i,"."))
#   i_set <- c()
#   iter_row <- allisolates_sets[i,]
#   for (j in seq_along(genes)){
#     gene <- genes[j]
#     val <- iter_row[[gene]]
#     if (val == 1){
#       i_set <- c(i_set, gene)
#     }
#   }
#   cat(paste0(i_set,"\n"))
#   if (length(i_set > 0)){
#     allisolates_sets$set[[i]] <- paste0(i_set,collapse = ",")
#   }
# }
# 
# allisolates_sets <- allisolates_sets %>% 
#   select(id, set) %>% 
#   mutate(set_vect = stringr::str_split(set, pattern=","),
#          mcr = allisolates$mcr,
#          match = FALSE,
#          match_lhs = FALSE)

allisolates_sets <- allisolates_sets %>%
  select(id, set, mcr) %>%
  mutate(match = FALSE,
         match_lhs = FALSE)

for (i in 1:nrow(allisolates_sets)){
  for (j in 1:nrow(rules)){
    # set <- unlist(allisolates_sets$set_vect[i])
    set <- unlist(allisolates_sets$set[i])
    rule_lhs <- unlist(rules$lhs[j])
    rule <- c(unlist(rules$lhs[j]), as.character(rules$rhs[j]))
    
    cat("Comparing set", i, "with rule", j, "\n")
    
    if(allisolates_sets$mcr[i] == 1){
      
      if(all(rule %in% set)){
        allisolates_sets$match[i] <- TRUE
        break
      } else {
        allisolates_sets$match[i] <- FALSE
      }
      
    } else {
      if(all(rule_lhs %in% set)){
        allisolates_sets$match_lhs[i] <- TRUE
        break
      } else {
        allisolates_sets$match_lhs[i] <- FALSE
      }
    }
    
  }
}

saveRDS(allisolates_sets, "allisolates_sets_comparison.RDS")

allisolates_sets_comparison <- readRDS("allisolates_sets_comparison.RDS") %>% 
  mutate(pctself = allisolates$pctself)
allisolates_sets_onlyARM <- readRDS("allisolates_sets_onlyARM.RDS") %>% 
  mutate(pctself = allisolates$pctself)
allisolates_sets_onlyLR <- readRDS("allisolates_sets_onlyLR.RDS") %>% 
  mutate(pctself = allisolates$pctself)
allisolates_sets_onlyGLM <- readRDS("allisolates_sets_onlyGLM.RDS") %>% 
  mutate(pctself = allisolates$pctself)
allisolates_sets_onlyADASYN <- readRDS("allisolates_sets_onlyADASYN.RDS") %>% 
  mutate(pctself = allisolates$pctself)

## Create match sets - All Rules
rule_matches_all <- allisolates_sets_comparison %>%
  filter(match == TRUE) %>% 
  mutate(comparison = "all rules",
         method = "ARM and LR",
         algorithm = "All",
         set = "rule matches")

rule_non_matches_all <- allisolates_sets_comparison %>%
  filter(match != TRUE) %>% 
  mutate(comparison = "all rules",
         method = "ARM and LR",
         algorithm = "All",
         set = "rule non-matches")

lhs_rule_matches_all <- allisolates_sets_comparison %>%
  filter(match_lhs == TRUE,
         mcr == 0) %>% 
  mutate(comparison = "all rules",
         method = "ARM and LR",
         algorithm = "All",
         set = "rule lhs matches")

lhs_rule_non_matches_all <- allisolates_sets_comparison %>%
  filter(match_lhs != TRUE,
         mcr == 0) %>% 
  mutate(comparison = "all rules",
         method = "ARM and LR",
         algorithm = "All",
         set = "rule lhs non-matches")

## Create match sets - ARM Rules
rule_matches_ARM <- allisolates_sets_onlyARM %>%
  filter(match == TRUE) %>% 
  mutate(comparison = "ARM rules",
         method = "ARM",
         algorithm = "apriori",
         set = "rule matches")

rule_non_matches_ARM <- allisolates_sets_onlyARM %>%
  filter(match != TRUE) %>% 
  mutate(comparison = "ARM rules",
         method = "ARM",
         algorithm = "apriori",
         set = "rule non-matches")

lhs_rule_matches_ARM <- allisolates_sets_onlyARM %>%
  filter(match_lhs == TRUE,
         mcr == 0) %>% 
  mutate(comparison = "ARM rules",
         method = "ARM",
         algorithm = "apriori",
         set = "rule lhs matches")

lhs_rule_non_matches_ARM <- allisolates_sets_onlyARM %>%
  filter(match_lhs != TRUE,
         mcr == 0) %>% 
  mutate(comparison = "ARM rules",
         method = "ARM",
         algorithm = "apriori",
         set = "rule lhs non-matches")

## Create match sets - LR Rules
rule_matches_LR <- allisolates_sets_onlyLR %>%
  filter(match == TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM and ADASYN",
         set = "rule matches")

rule_non_matches_LR <- allisolates_sets_onlyLR %>%
  filter(match != TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM and ADASYN",
         set = "rule non-matches")

lhs_rule_matches_LR <- allisolates_sets_onlyLR %>%
  filter(match_lhs == TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM and ADASYN",
         set = "rule lhs matches")

lhs_rule_non_matches_LR <- allisolates_sets_onlyLR %>%
  filter(match_lhs != TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM and ADASYN",
         set = "rule lhs non-matches")

## Create match sets - GLM Rules
rule_matches_GLM <- allisolates_sets_onlyGLM %>%
  filter(match == TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM",
         set = "rule matches")

rule_non_matches_GLM <- allisolates_sets_onlyGLM %>%
  filter(match != TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM",
         set = "rule non-matches")

lhs_rule_matches_GLM <- allisolates_sets_onlyGLM %>%
  filter(match_lhs == TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM",
         set = "rule lhs matches")

lhs_rule_non_matches_GLM <- allisolates_sets_onlyGLM %>%
  filter(match_lhs != TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "GLM",
         set = "rule lhs non-matches")

## Create match sets - ADASYN Rules
rule_matches_ADASYN <- allisolates_sets_onlyADASYN %>%
  filter(match == TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "ADASYN",
         set = "rule matches")

rule_non_matches_ADASYN <- allisolates_sets_onlyADASYN %>%
  filter(match != TRUE) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "ADASYN",
         set = "rule non-matches")

lhs_rule_matches_ADASYN <- allisolates_sets_onlyADASYN %>%
  filter(match_lhs == TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "ADASYN",
         set = "rule lhs matches")

lhs_rule_non_matches_ADASYN <- allisolates_sets_onlyADASYN %>%
  filter(match_lhs != TRUE,
         mcr == 0) %>% 
  mutate(comparison = "LR rules",
         method = "LR",
         algorithm = "ADASYN",
         set = "rule lhs non-matches")


rule_matches_vs_rule_non_matches_plotdata <- rbind(rule_matches_all,
                                                   rule_non_matches_all,
                                                   lhs_rule_matches_all,
                                                   lhs_rule_non_matches_all,
                                                   rule_matches_ARM,
                                                   rule_non_matches_ARM,
                                                   lhs_rule_matches_ARM,
                                                   lhs_rule_non_matches_ARM,
                                                   rule_matches_LR,
                                                   rule_non_matches_LR,
                                                   lhs_rule_matches_LR,
                                                   lhs_rule_non_matches_LR,
                                                   rule_matches_GLM,
                                                   rule_non_matches_GLM,
                                                   lhs_rule_matches_GLM,
                                                   lhs_rule_non_matches_GLM,
                                                   rule_matches_ADASYN,
                                                   rule_non_matches_ADASYN,
                                                   lhs_rule_matches_ADASYN,
                                                   lhs_rule_non_matches_ADASYN)
saveRDS(rule_matches_vs_rule_non_matches_plotdata, "rule_matches_vs_rule_non_matches_plotdata.RDS")
write.csv(rule_matches_vs_rule_non_matches_plotdata %>% mutate(set_vect = paste0(set_vect)),
          file = "rule_matches_vs_rule_non_matches_plotdata.csv")

### SET COMPARISONS - STATS TESTS (ALL)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "greater")

## lhs Rule matches vs. non-matches (without MCR)
# wilcox.test(pctself ~ set,
#             data = rule_matches_vs_rule_non_matches_plotdata %>% 
#               filter(set %in% c("rule lhs matches", "rule lhs non-matches"),
#                      algorithm == "All"),
#             paired = FALSE)

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "All"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (ARM)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "apriori"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (GLM)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "GLM"),
            paired = FALSE,
            alternative = "less")

### SET COMPARISONS - STATS TESTS (ADASYN)
## Rule matches vs. non-matches
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "greater")

## Rule matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule matches", "rule lhs matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "less")

## Rule non-matches: mcr vs. non-mcr (lhs)
wilcox.test(pctself ~ set,
            data = rule_matches_vs_rule_non_matches_plotdata %>% 
              filter(set %in% c("rule non-matches", "rule lhs non-matches"),
                     algorithm == "ADASYN"),
            paired = FALSE,
            alternative = "less")

################
## Comparisons for Gene Sets in the Rules (from ARM and GLM)
rules_sets <- readRDS("../web_crawler/all_rules.RDS") %>% 
  mutate(pct_self = NA)


for (i in 1:nrow(rules_sets)){
  pct_self <- length(which(unlist(rules_sets$rule[i]) %in% self_vector))/(rules_sets$size[i])
  rules_sets$pct_self[i] <- pct_self
}

saveRDS(rules_sets, "../web_crawler/all_rules.RDS")
readr::write_csv(rules_sets %>% select(-rule), "../web_crawler/all_rules.csv")
