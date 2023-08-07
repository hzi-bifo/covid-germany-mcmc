# R script to adjust the metadata into a format that can be used as input for the Lineage SDplots

library("stringr")
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#print(args[1])
#setwd("/home/sreimering/Documents/Coronavirus/Hadi/importation_lineages_DE_casesampling/data")
setwd(args[1])


file <- "clusterSamples_DTA_MCC_0.5.tsv"
output_file <- "clusterSamples_DTA_MCC_0.5_reformatted.tsv"
data <- read.table(file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#todo: rename date column, rename lineage column, add country column, add state column, add submitting and originating lab
# rename columns
names(data)[names(data) == 'Collection.date'] <- 'date'
#names(data)[names(data) == 'cluster'] <- 'pangolin_lineage'
names(data)[names(data) == 'Accession.ID'] <- 'gisaid_epi_isl'
names(data)[names(data) == 'Host'] <- 'host'

# add new columns
data$pangolin_lineage <- paste(sapply(strsplit(data$cluster, '-'), function(x) x[3]), sapply(strsplit(data$cluster, '_'), function(x) x[4]), sep = "_")
data$country <- "Germany"
data$division <- sapply(strsplit(data$Location, '/'), function(x) str_trim(x[3]))
data$originating_lab <- NA
data$submitting_lab <- NA

# rename states
data$division[data$division == "Baden-Wurttemberg"] <- "Baden-Wuerttemberg"
data$division[data$division == "Baden-Württemberg"] <- "Baden-Wuerttemberg"
data$division[data$division == "Lower Saxony"] <- "Niedersachsen"
data$division[data$division == "Mecklenburg-Western Pomerania"] <- "Mecklenburg-Vorpommern"
data$division[data$division == "North Rhine-Westphalia"] <- "North Rhine Westphalia"
data$division[data$division == "Rhineland-Palatinate"] <- "Rheinland-Pfalz"

write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
