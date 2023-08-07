#Rscript
library(stringr)
args = commandArgs(trailingOnly=TRUE)
state <- args[1]
metadata_file <- args[2]
out_list_file <- args[3]

metadata <- read.table(metadata_file, sep="\t", head=TRUE, na.strings=c("NA", ""), fill=TRUE, stringsAsFactors=FALSE, quote="|")
metadata <- metadata[metadata$Virus.name != "Virus.name",]
metadata$country <- sapply(strsplit(metadata$Location, '/'), function(x) str_trim(x[2]))
write.table(metadata$Accession.ID[metadata$country == state], out_list_file, row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")
