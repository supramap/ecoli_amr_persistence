############################################################
## Script for generating Logisitc Regression              ##
## for E. coli genotype sets containing `mcr`             ##
## By: Gabriel Zenarosa, Ph.D., Colby T. Ford, Ph.D.,     ##
##     David Brown, Kevin Smith, and Daniel Janies, Ph.D. ##
############################################################

library(jsonlite)
library(stringr)
library(glmnet)
#library(caret)
#library(arules)

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
gdf.mcr$mcr <- as.factor(gdf.mcr$mcr)
#for(col in 1:ncol(gdf.mcr)){gdf.mcr[,col] <- as.factor(gdf.mcr[,col])}
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

#### GLMNET
## Setup Parallel Backend
library(doMC)
n_cores <- parallel::detectCores()-1
cat("Number of cores:", n_cores)
registerDoMC(n_cores)

x <- Matrix(model.matrix(f, gdf.mcr)[, -1], sparse = T)
y <- Matrix(gdf.mcr$mcr, sparse = T)
set.seed(123456789)
model.cv <- cv.glmnet(x = x,
                      y = y,
                      family = "binomial",
                      alpha = 1,
                      #nfolds = floor(length(y) / 10),
                      type.measure = "auc",
                      parallel = T)
saveRDS(model.cv, "mcr_lr_model.RDS")
model.cv <- readRDS("mcr_lr_model.RDS")

#library(DMwR)
#gdf.mcr.smote <- SMOTE(mcr ~ .,
#                       gdf.mcr,
#                       perc.over = 100 * ceiling(nrow(gdf.mcr[gdf.mcr$mcr == 0,]) / nrow(gdf.mcr[gdf.mcr$mcr == 1,]) - 1),
#                       perc.under = 100 + (100 * nrow(gdf.mcr[gdf.mcr$mcr == 1,]) / nrow(gdf.mcr[gdf.mcr$mcr == 0,])))
library(smotefamily)
gdf.mcr.adasyn <- ADAS(gdf.mcr[,-ncol(gdf.mcr)], gdf.mcr$mcr)
gdf.mcr.smote <- gdf.mcr.adasyn$data
colnames(gdf.mcr.smote)[colnames(gdf.mcr.smote) == "class"] <- "mcr"
#table(gdf.mcr$mcr)
#table(gdf.mcr.smote$mcr)
x.smote <- Matrix(model.matrix(f, gdf.mcr.smote)[, -1], sparse = T)
y.smote <- Matrix(gdf.mcr.smote$mcr, sparse = T)
set.seed(123456789)
model.cv.adasyn <- cv.glmnet(x = x.smote,
                            y = y.smote,
                            family = "binomial",
                            alpha = 1,
                            nfolds = floor(length(y) / 10),
                            type.measure = "auc",
                            parallel = T)
saveRDS(model.cv.adasyn, "mcr_lr_model_adasyn.RDS")
model.cv.adasyn <- readRDS("mcr_lr_model_adasyn.RDS")

### Use adasyn/smote?
#model.cv <- model.cv.smote

coef.mcr <- coef(model.cv, s = "lambda.min")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
#write.csv(as.data.frame(exp(cbind(Odds.Ratio=coef.mcr.nz))), "oddsRatioInt.csv")
View(data.frame(AMR=coef.mcr.names[nonzeroCoef(coef.mcr)], Odds.Ratio=exp(coef.mcr[nonzeroCoef(coef.mcr)])))

## Look at all prevalidated fits to determine best Lambda, based on Sensitivity
# allfits <- as.data.frame(model.cv$fit.preval[,1:length(model.cv$lambda)])
# colnames(allfits) <- paste0("Lambda:",model.cv$lambda)
# allfits$actual <- gdf.mcr$mcr

# cor.test(as.numeric(gdf.mcr$mcr), allfits[,1])
# MLmetrics::Sensitivity(y_pred = allfits[,1],
#                        y_true = gdf.mcr$mcr,
#                        positive = lev[1])

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
# New data From Jan 9, 2019 to May 8, 2019
# json.stream <- fromJSON("https://www.ncbi.nlm.nih.gov/pathogens/ngram?start=0&limit=1000000&q=%5Bdisplay()%2Chist(geo_loc_name%2Cisolation_source%2Cepi_type%2Ccollected_by%2Chost%2Cproperty%2Ctarget_creation_date)%5D.from(pathogen).usingschema(%2Fschema%2Fpathogen).matching(kmer_group%3D%3D%5B%22PDG000000004.1238%22%5D+and+q%3D%3D%22taxgroup_name%253A%2522E.coli%2520and%2520Shigella%2522%22).sort(target_creation_date%2Cdesc)&_search=false&rows=20&page=1&sidx=target_creation_date&sord=desc")
# saveRDS(json.stream, file = "e.coli_validation.RDS")
newjson.stream <- readRDS("e.coli_validation.RDS")

newid <- substr(newjson.stream[["ngout"]][["data"]][["content"]][["id"]], 19, 33)
newcreation_date_time <- as.POSIXct(newjson.stream[["ngout"]][["data"]][["content"]][["target_creation_date"]], format = "%Y-%m-%dT%H:%M:%SZ")
newcollection_year <- as.numeric(substr(newjson.stream[["ngout"]][["data"]][["content"]][["collection_date"]], 1, 4))
newcollection_year[is.na(newcollection_year)] <- as.numeric(format(newcreation_date_time[is.na(newcollection_year)], "%Y"))
newlocation <- newjson.stream[["ngout"]][["data"]][["content"]][["geo_loc_name"]]
newisolation_type <- newjson.stream[["ngout"]][["data"]][["content"]][["epi_type"]]
newisolation_source <- newjson.stream[["ngout"]][["data"]][["content"]][["isolation_source"]]
newdf <- data.frame(newid, newcreation_date_time, newcollection_year, newlocation, newisolation_type, newisolation_source)
newgenotypes <- newjson.stream[["ngout"]][["data"]][["content"]][["AMR_genotypes"]]
names(newgenotypes) <- newid

newgenotypes <- newgenotypes[newcreation_date_time >= as.POSIXct("2019-01-09 00:00:00", format = "%Y-%m-%d %H:%M:%S")]
newgenoid <- newid[newcreation_date_time >= as.POSIXct("2019-01-09 00:00:00", format = "%Y-%m-%d %H:%M:%S")]


# ---------
# Genotypes
# ---------
#genotypesall <- unique(sort(unlist(genotypes)))
newgenobinv <- lapply(newgenotypes, function(x){(genotypesall %in% x) * 1})
newgdf <- as.data.frame(matrix(unlist(newgenobinv), ncol = length(genotypesall), byrow = T))
names(newgdf) <- genotypesall
newgdf.mcr <- newgdf[,-grep("mcr", names(newgdf))]
newgdf.mcr$mcr <- (rowSums(newgdf[,grep("mcr", names(newgdf))]) > 0) * 1
newgdf.mcr$mcr <- as.factor(newgdf.mcr$mcr)
# for(col in 1:ncol(newgdf.mcr)){newgdf.mcr[,col] <- as.factor(newgdf.mcr[,col])}

# -------------
# Genotype Sets
# -------------
# genotypes.mcr <- genotypes
# genotypes.mcr.idx <- grep("mcr", genotypes)
# genotypes.mcr[genotypes.mcr.idx] <- lapply(genotypes[genotypes.mcr.idx], function(x){x[-grep("mcr", x)]})
# genotypeset.mcr <- lapply(genotypes.mcr,
#                           function(x){
#                             # add ticks (`) only as needed
#                             y <- format(as.formula(paste0(ifelse(length(x) > 0, paste0("`", paste(x, collapse = "`:`"), "`"), ""), "~1")))[1]
#                             return(substr(y, 1, regexpr("~", y)[1] - 2))
#                           }
# )
# genotypesets.mcr <- unique(genotypeset.mcr)
#genotypesets.mcr.gt1 <- genotypesets.mcr[grep(":", genotypesets.mcr)]
#f <- as.formula(paste("mcr ~ 1+.", paste(genotypesets.mcr.gt1, collapse = "+"), sep = "+"))
newx <- Matrix(model.matrix(f, newgdf.mcr)[, -1], sparse = T)
newy <- Matrix(newgdf.mcr$mcr, sparse = T)

## Run validation
test.sv <- predict(model.cv, newx = x, s = "lambda.min", type = "response")
test.cv <- predict(model.cv, newx = newx, s = "lambda.min", type = "response")

test.sv.adasyn <- predict(model.cv.adasyn, newx = x, s = "lambda.min", type = "response")
test.cv.adasyn <- predict(model.cv.adasyn, newx = newx, s = "lambda.min", type = "response")


# table((test.sv > 0.15) * 1, gdf.mcr$mcr)
# table((test.cv > 0.01) * 1, newgdf.mcr$mcr)

######################
## Get Odds Ratios for Each Model
model.cv <- readRDS("mcr_lr_model.RDS")
coef.mcr <- coef(model.cv, s = "lambda.min")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
oddratios_glm <- data.frame(AMR=coef.mcr.names[nonzeroCoef(coef.mcr)], Odds.Ratio=exp(coef.mcr[nonzeroCoef(coef.mcr)]))

model.cv.adasyn <- readRDS("mcr_lr_model_adasyn.RDS")
model.cv <- model.cv.adasyn
coef.mcr <- coef(model.cv, s = "lambda.min")
coef.mcr.names <- rownames(coef.mcr)
coef.mcr.nz <- Matrix(coef.mcr[nonzeroCoef(coef.mcr)], sparse = T)
rownames(coef.mcr.nz) <- coef.mcr.names[nonzeroCoef(coef.mcr)]
oddratios_adasyn <- data.frame(AMR=coef.mcr.names[nonzeroCoef(coef.mcr)], Odds.Ratio=exp(coef.mcr[nonzeroCoef(coef.mcr)]))

oddsratios <- full_join(oddratios_glm, oddratios_adasyn, by = "AMR")
colnames(oddsratios) <- c("GeneSet", "OddsRatio_glm", "OddsRatio_adasyn")
write.csv(oddsratios, "mcr_lr_model_oddsratios.csv", row.names = FALSE)

######################
## Generate ROC curves for Each Model
library(pROC)
library(ggplot2)

model.cv.glm <- readRDS("mcr_lr_model.RDS")
model.cv.adasyn <- readRDS("mcr_lr_model_adasyn.RDS")

test.cv.glm <- predict(model.cv.glm, newx = newx, s = "lambda.min", type = "response")
test.cv.adasyn <- predict(model.cv.adasyn, newx = newx, s = "lambda.min", type = "response")

rocobj.glm <- roc(newgdf.mcr$mcr, test.cv.glm[,1])
rocobj.adasyn <- roc(newgdf.mcr$mcr, test.cv.adasyn[,1])

ggroc(list(glm=rocobj.glm, adasyn=rocobj.adasyn)) +
  theme_minimal() +
  ggtitle("ROC curve") + 
  xlab("Specificity") +
  ylab("Sensitivity") + 
  labs(fill = "Model") + 
  geom_segment(aes(x = 1,
                   xend = 0,
                   y = 0,
                   yend = 1),
               color="grey",
               linetype="dashed")

## Analyze and compares the ROC curves
auc(rocobj.glm)
auc(rocobj.adasyn)

roc.test(rocobj.glm, rocobj.adasyn)

allcoords_glm <- coords(rocobj.glm,
                        ret = "all",
                        transpose = FALSE) %>% 
  mutate(is_best = ifelse(.$threshold == coords(rocobj.glm, "best", ret = c("threshold"), transpose = FALSE), TRUE, FALSE),
         model = "glm")

allcoords_adasyn <- coords(rocobj.adasyn,
                          ret = "all",
                          transpose = FALSE) %>% 
  mutate(is_best = ifelse(.$threshold == coords(rocobj.adasyn, "best", ret = c("threshold"), transpose = FALSE), TRUE, FALSE),
         model = "adasyn")

allcoords <- rbind(allcoords_glm, allcoords_adasyn)

write.csv(allcoords, "ROCcoordinates.csv", row.names = FALSE)

## Create Precision-Recall Curve
library(DMwR)

## GLM Model
DMwR::PRcurve(test.sv, gdf.mcr$mcr)
DMwR::PRcurve(test.cv, newgdf.mcr$mcr)

## Adasyn
DMwR::PRcurve(test.sv, gdf.mcr$mcr)
DMwR::PRcurve(test.cv.adasyn, newgdf.mcr$mcr)
