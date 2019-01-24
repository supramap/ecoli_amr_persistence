library(jsonlite)
library(stringr)
library(glmnet)

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
genotypes <- fgenotypes

# ---------
# Genotypes
# ---------
genotypesall <- unique(unlist(genotypes))
genobinv <- lapply(genotypes, function(x){(genotypesall %in% x) * 1})
gdf <- as.data.frame(matrix(unlist(genobinv), ncol = length(genotypesall), byrow = T))
names(gdf) <- genotypesall
gdf.mcr <- gdf[,-grep("mcr", names(gdf))]
gdf.mcr$mcr <- (rowSums(gdf[,grep("mcr", names(gdf))]) > 0) * 1
#model.null <- glm(mcr ~ 1, data = gdf.mcr, family = binomial)
#model.full <- glm(mcr ~ ., data = gdf.mcr, family = binomial)
#mixed <- step(model.null,
#              scope = list(lower = formula(model.null), upper = formula(model.full)),
#              direction = "forward")
#x <- Matrix(as.matrix(gdf.mcr[,-which(names(gdf.mcr) == "mcr")]), sparse = T)
f <- as.formula("mcr ~ .+1")
x <- Matrix(model.matrix(f, gdf.mcr)[, -1], sparse = T)
y <- Matrix(gdf.mcr$mcr, sparse = T)
#model.full <- glmnet(x = x, y = y, family = "binomial", lambda = 0)
model.cv <- cv.glmnet(x = x, y = y, family = "binomial", type.measure = "class")
#plot(model.cv)
#model.full <- glmnet(x = x, y = y, family = "binomial", lambda = model.cv$lambda.1se)
#model.full <- model.cv$glmnet.fit
#exp(cbind(Odds.Ratio=coef(model.full), Confidence.Interval=confint(model.full)))
coef.mcr <- coef(model.cv, s = "lambda.1se")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
#logistic.display(model.full)

coef.all <- unlist(lapply(genotypesall, function(x){coef(glm(get(x) ~ 1, family = "binomial", data = gdf))}))
names(coef.all) <- genotypesall

# -------------
# Genotype Sets
# -------------
genotypes.mcr <- genotypes
genotypes.mcr.idx <- grep("mcr", genotypes.mcr)
genotypes.mcr[genotypes.mcr.idx] <- lapply(genotypes[genotypes.mcr.idx], function(x){x[-grep("mcr", x)]})
genotypeset.mcr <- lapply(genotypes.mcr, paste, collapse = "`:`")
genotypesets.mcr <- paste0("`", unique(unlist(genotypeset.mcr)), "`")
genotypesets.mcr <- genotypesets.mcr[-grep("``", genotypesets.mcr)]
f <- as.formula(paste("mcr ~ .+1", paste(genotypesets.mcr, collapse = "+"), sep = "+"))
x <- Matrix(model.matrix(f, gdf.mcr)[, -1], sparse = T)
y <- Matrix(gdf.mcr$mcr, sparse = T)
model.cv <- cv.glmnet(x = x, y = y, family = "binomial", type.measure = "class")
coef.mcr <- coef(model.cv, s = "lambda.1se")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
exp(cbind(Odds.Ratio=coef.mcr.nz))
#write.csv(as.data.frame(exp(cbind(Odds.Ratio=coef.mcr.nz))), "oddsRatioInt.csv")
