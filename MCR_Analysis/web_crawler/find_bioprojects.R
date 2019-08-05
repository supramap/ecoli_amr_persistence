#######################
## Script for finding all BioProject numbers for publication crawling
## By: Colby T. Ford, Ph.D., Gabriel Zenarosa, Ph.D., David Brown, Kevin Smith, and Daniel Janies, Ph.D.
#######################

## Load in Packages
library(jsonlite)
library(tidyr)
library(dplyr)

## Generate New Extract from the NCBI Database
#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=100000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Cdesc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=desc)")
#saveRDS(json.stream, file = "json.stream.RDS")

## Read in Data
json.stream <- readRDS("e.coli.RDS")

## Get IDs
id <- json.stream[["ngout"]][["data"]][["content"]][["id"]] %>% as.data.frame()
colnames(id) <- "id"

## Separate IDs into PDG and PDT numbers
ids <- id %>% separate(id, c("pdg", "pdt"), sep = "_")

## Get BioProject numbers
bioprojects <- json.stream[["ngout"]][["data"]][["content"]][["bioproject_acc"]] %>% as.data.frame()
colnames(bioprojects) <- "bioproject"

bioprojects <- data.frame(ids, bioprojects)
bioprojects$bioproject_url <- paste0("https://www.ncbi.nlm.nih.gov/bioproject/",bioprojects$bioproject)

## Join to MCR Groups (from Odds Ratio Analysis)
# mcr_groups <- read.csv("MCR_Groups_PDTs.csv")
mcr_groups <- read.csv("allPDTs.csv")

bioprojects <- inner_join(mcr_groups, bioprojects)

# write.csv(bioprojects, "BioProjects_2019-14-02.csv")

######################
## Find all citations
library(rvest)
library(stringr)
library(reshape)

bioprojects_unique <- unique(bioprojects$bioproject)
publications <- data.frame(bioproject = character(0),
                           publications = character(0))

for (i in 1:length(bioprojects_unique)){
  bioproject <- bioprojects_unique[i]
  # content <- read_html(paste0("https://www.ncbi.nlm.nih.gov/bioproject/",
  #                             bioprojects_unique[i])) %>% 
  #   html_nodes('#CombinedTable') %>% 
  #   html_table() %>% 
  #   .[[1]] %>% 
  #   filter(X1 == "Publications") %>% 
  #   select(X2)
  
  content <- read_html(paste0("https://www.ncbi.nlm.nih.gov/bioproject/",
                              bioprojects_unique[i])) %>% 
    html_nodes("a.RegularLink") %>% 
    html_attr("href") %>% 
    as.data.frame() %>% 
    filter(str_detect(., "/pubmed/")) %>% 
    .[[1]] %>% 
    as.data.frame()
  
  colnames(content) <- "publications"
  
  if (length(content$publications) > 0){
    iteroutput <- data.frame(bioproject = bioproject,
                             publications = content)
    
    publications <- rbind(publications, iteroutput)
  }

  cat(paste0(i, ". BioProject: ", as.character(bioprojects_unique[i]), "\n"))
}

#publications_orig <- publications
## Derive URLS for Publications
publications$publication_url <- paste0("https://www.ncbi.nlm.nih.gov", publications$publications)

## Join Publications to PDT numbers via BioProject
publication_output <- inner_join(publications, bioprojects)

publications_by_mcr_group <- publication_output %>%
  select(-c("pdt", "pdg", "publications", "bioproject")) %>% 
  unique() %>% 
  cast(publication_url~mcr_group)

# write.csv(publications_by_mcr_group, "publications_by_mcr_group.csv")
write.csv(publications_by_mcr_group, "publications_by_ARM_mcr_group.csv")
