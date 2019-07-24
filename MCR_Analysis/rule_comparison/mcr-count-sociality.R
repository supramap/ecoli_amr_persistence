library(tidyr)
library(dplyr)

rules <- readRDS("mcr_rules_all.RDS")
gene <- read.csv("AMR_FunctionalMechanisms.csv", header=T)

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

all_rules <- readRDS("mcr_rules_all.RDS")
val_rules <- readRDS("mcr_validation_rules.RDS")

#Join excluding column-list called "rule" (column 8) because of error and irrelevancy
inner_rules <- inner_join(all_rules[,-8], val_rules[,-8], by="lhs")
left_rules <- left_join(all_rules[,-8], val_rules[,-8], by="lhs")

#Set non-matching sets conf,supp,lift to 0
left_rules$support.y[is.na(left_rules$support.y)] <- 0
left_rules$confidence.y[is.na(left_rules$confidence.y)] <- 0
left_rules$lift.y[is.na(left_rules$lift.y)] <- 0

#RMSD Calculations
rmsd_inner_support <- sqrt((sum((abs(inner_rules$support.x - inner_rules$support.y))^2))/nrow(inner_rules))
rmsd_left_support <- sqrt((sum((abs(left_rules$support.x - left_rules$support.y))^2))/nrow(left_rules))

rmsd_inner_confidence <- sqrt((sum((abs(inner_rules$confidence.x - inner_rules$confidence.y))^2))/nrow(inner_rules))
rmsd_left_confidence <- sqrt((sum((abs(left_rules$confidence.x - left_rules$confidence.y))^2))/nrow(left_rules))

rmsd_inner_lift <- sqrt((sum((abs(inner_rules$lift.x - inner_rules$lift.y))^2))/nrow(inner_rules))
rmsd_left_lift <- sqrt((sum((abs(left_rules$lift.x - left_rules$lift.y))^2))/nrow(left_rules))
