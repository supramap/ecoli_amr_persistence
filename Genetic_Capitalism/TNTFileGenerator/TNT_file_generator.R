library(jsonlite)
library(stringr)

# Run these *nix commands
#sed 's/'\''//g' PDG000000004.1024.reference_target.tree.newick | sed 's/:-*[0-9]\.*[0-9]*\(e-[0-9]*\)*//g' | sed 's/,'\('/'\('/g' | sed 's/,/ /g' > e.coli.paren
#sed 's/'\)'*'\('*'\ '*PDT/,PDT/g' e.coli.paren | sed 's/^,//g' | sed 's/'\)'*;$//g' > e.coli.ids
#---------
pdt <- unique(names(read.csv("e.coli.ids")))

#json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=1000000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Casc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=asc)")
#saveRDS(json.stream, file = "e.coli.RDS")
json.stream <- readRDS("e.coli.RDS")

id <- substr(json.stream[["ngout"]][["data"]][["content"]][["id"]], 19, 33)
creation_date_time <- as.POSIXct(json.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")
collection_year <- as.numeric(substr(json.stream[["ngout"]][["data"]][["content"]][["collection_date"]], 1, 4))
collection_year[is.na(collection_year)] <- as.numeric(format(creation_date_time[is.na(collection_year)], "%Y"))
location <- json.stream[["ngout"]][["data"]][["content"]][["geo_loc_name"]]
isolation_type <- json.stream[["ngout"]][["data"]][["content"]][["epi_type"]]
isolation_source <- json.stream[["ngout"]][["data"]][["content"]][["isolation_source"]]
genotypes <- json.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
df <- data.frame(id, creation_date_time, collection_year, location, isolation_type, isolation_source)
df <- df[df$id %in% pdt,]
names(genotypes) <- id
genotypes[!(names(genotypes) %in% pdt)] <- NULL
#genotypes[genotypes == "NULL"] <- NULL

# ---------
# Genotypes
# ---------
genotypesall <- unique(unlist(genotypes))
genobinv <- lapply(genotypes, function(x){(genotypesall %in% x) * 1})
genotype_count <- unlist(lapply(genobinv, sum))
genobins <- lapply(genobinv, paste, collapse = "")
dfbins <- data.frame(df, genotype_count, genobins=unlist(genobins))
dfbins <- dfbins[order(dfbins$genotype_count, dfbins$collection_year, dfbins$creation_date_time),]
genobins <- dfbins$genobins
names(genobins) <- dfbins$id

cat("mxram 800000", "xread", paste(c(length(genotypesall), length(genobins)), collapse = " "), file = "e.coli.genotypes.tnt.head", sep = "\n")
cat(as.vector(paste(names(genobins), unlist(genobins), sep = " ")), file = "e.coli.genotypes.tnt.head", sep = "\n", append = T)
cat(";", "tread", file = "e.coli.genotypes.tnt.head", sep = "\n", append = T)

cat(";", "", "cnames", file = "e.coli.genotypes.tnt.tail", sep = "\n")
for (i in 1:length(genotypesall)) {
  cat(paste0("{", (i - 1), " 'genotype_", genotypesall[i], "' absent present;"), file = "e.coli.genotypes.tnt.tail", sep = "\n", append = T)
}
cat(";", "proc /;", file = "e.coli.genotypes.tnt.tail", sep = "\n", append = T)

system("cat e.coli.genotypes.tnt.head e.coli.paren e.coli.genotypes.tnt.tail > e.coli.genotypes.tnt")
write.csv(dfbins, "e.coli.genotypes.meta.csv")

# ----------------------------
# Genotype Odds and Odds Ratio
# ----------------------------
dfbinv <- data.frame(matrix(unlist(genobinv), nrow=length(genobinv), byrow=T))
colnames(dfbinv) <- genotypesall
dfbinv <- data.frame(amr=rep(1, nrow(dfbinv)), dfbinv)
library(glmnet)
model <- cv.glmnet(x = as.matrix(dfbinv), y = dfbinv$amr, family = "binomial", type.measure = "class")

colSums(dfbinv)['blaEC']
genotypeodds <- rep(0.0, length(genotypesall))
names(genotypeodds) <- genotypesall
for (g in genotypesall) {
  genotypeodds[g] <- coef(glm(get(g) ~ blaEC, family = binomial, data = dfbinv))['blaEC']
}
genotypeodds['blaEC'] <- coef(glm(blaEC ~ 1, family = binomial, data = dfbinv))
write.csv(as.data.frame(exp(cbind(Odds.Ratio=genotypeodds))), "oddsRatioSingles.csv")
# dfbinv0 <- dfbinv[dfbinv$blaEC == 0,]
# genotypeodds0 <- rep(0.0, length(genotypesall))
# names(genotypeodds0) <- genotypesall
# for (g in genotypesall) {
#   genotypeodds0[g] <- coef(glm(get(g) ~ 1, family = binomial, data = dfbinv0))
# }

# -------------
# Genotype Sets
# -------------
genotypeset <- lapply(genotypes, paste, collapse="_")
genotypesets <- unique(unlist(genotypeset))
genobinsetv <- lapply(genotypeset, function(x){(genotypesets == x) * 1})
genotypeset_count <- unlist(lapply(genobinsetv, sum))
genobinset <- lapply(genobinsetv, paste, collapse = "")
dfbinset <- data.frame(df, genotypeset_count, genobinset=unlist(genobinset))
dfbinset <- dfbinset[order(dfbinset$genotypeset_count, dfbinset$collection_year, dfbinset$creation_date_time),]
genobinset <- dfbinset$genobinset
names(genobinset) <- dfbinset$id

cat("mxram 800000", "xread", paste(c(length(genotypesets), length(genobinset)), collapse = " "), file = "e.coli.tnt.head", sep = "\n")
cat(as.vector(paste(names(genobinset), unlist(genobinset), sep = " ")), file = "e.coli.tnt.head", sep = "\n", append = T)
cat(";", "tread", file = "e.coli.tnt.head", sep = "\n", append = T)

cat(";", "", "cnames", file = "e.coli.tnt.tail", sep = "\n")
for (i in 1:length(genotypesets)) {
  cat(paste0("{", (i - 1), " 'genotype_set_", genotypesets[i], "' absent present;"), file = "e.coli.tnt.tail", sep = "\n", append = T)
}
cat(";", "proc /;", file = "e.coli.tnt.tail", sep = "\n", append = T)

system("cat e.coli.tnt.head e.coli.paren e.coli.tnt.tail > e.coli.tnt")
write.csv(dfbinset, "e.coli.meta.csv")
