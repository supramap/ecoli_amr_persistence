#######################
## Script for generating Logisitc Regression for E. coli genotype sets containing `mcr`
## By:  Gabriel Zenarosa, Ph.D., Colby T. Ford, Ph.D., David Brown, Kevin Smith, and Daniel Janies, Ph.D.
#######################

library(jsonlite)
library(stringr)
library(glmnet)
library(arules)

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
# f <- as.formula("mcr ~ 1+.")
# x <- Matrix(model.matrix(f, gdf.mcr)[, -1], sparse = T)
# y <- Matrix(gdf.mcr$mcr, sparse = T)
# model.cv <- cv.glmnet(x = x, y = y, family = "binomial", type.measure = "class")
# #plot(model.cv)
# coef.mcr <- coef(model.cv, s = "lambda.1se")
# coef.mcr.names <- rownames(coef.mcr)
# coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
# rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]

# -------------
# Genotype Sets
# -------------
genotypes.mcr <- genotypes
genotypes.mcr.idx <- grep("mcr", genotypes)
genotypes.mcr[genotypes.mcr.idx] <- lapply(genotypes[genotypes.mcr.idx], function(x){x[-grep("mcr", x)]})
genotypeset.mcr <- lapply(genotypes.mcr,
  function(x){
    # add ticks (`) only as needed
    y <- format(as.formula(paste0(ifelse(length(x) > 0, paste0("`", paste(x, collapse = "`:`"), "`"), ""), "~1")))[1]
    return(substr(y, 1, regexpr("~", y)[1] - 2))
  }
)
genotypesets.mcr <- unique(genotypeset.mcr)
genotypesets.mcr.gt1 <- genotypesets.mcr[grep(":", genotypesets.mcr)]
f <- as.formula(paste("mcr ~ 1+.", paste(genotypesets.mcr.gt1, collapse = "+"), sep = "+"))
x <- Matrix(model.matrix(f, gdf.mcr)[, -1], sparse = T)
y <- Matrix(gdf.mcr$mcr, sparse = T)
model.cv <- cv.glmnet(x = x, y = y, family = "binomial", type.measure = "class", alpha = 0)
coef.mcr <- coef(model.cv, s = "lambda.1se")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
#write.csv(as.data.frame(exp(cbind(Odds.Ratio=coef.mcr.nz))), "oddsRatioInt.csv")

# ------------------
# Genotype Supersets
# ------------------
pdt.mcr.nz <- lapply(rownames(coef.mcr.nz)[-1], function(y){genoid[rowMeans(sapply(c("mcr", gsub('`', '', strsplit(y, ":")[[1]])), function(x){gdf.mcr[[x]]})) == 1]})
names(pdt.mcr.nz) <- rownames(coef.mcr.nz)[-1]
coef.mcr.nz.df <- data.frame(Odds.Ratio=exp(coef.mcr.nz)[-1, 1], Isolates=lengths(pdt.mcr.nz), URL=paste0("https://www.ncbi.nlm.nih.gov/pathogens/isolates/#/search/target_acc:", sapply(pdt.mcr.nz, paste, collapse = "%20OR%20target_acc:")))
# write.csv(coef.mcr.nz.df, "oddsRatio-MCR-vs-Sets.csv")

# ------------------
# Model Validation
# ------------------
# New data From Jan 9 to today
json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=1000000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(status%3D%3D%5B%22current%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Casc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=asc)")

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
genoid <- id[creation_date_time >= as.POSIXct("2019-09-01 20:14:38", format = "%Y-%m-%d %H:%M:%S")]


#Run validation
predict(model.cv, newx = [??] , s = "lambda.1se")