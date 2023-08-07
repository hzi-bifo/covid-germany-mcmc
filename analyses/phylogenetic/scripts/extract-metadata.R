#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#samples_file <- 'results/gisaid-20210428-samples.txt'
#metadata_file <- '../../data/phylogenetic/gisaid-20210503-metadata.tsv'
#metadata_sampled_file <- 'results/gisaid-20210428-metadata-sampled.txt'
samples_file <- args[1]
metadata_file <- args[2]
metadata_sampled_file <- args[3]

metadata <- read.table(metadata_file, sep="\t", head=TRUE, na.strings=c("NA", ""), fill=TRUE, stringsAsFactors=FALSE, quote="|")
samples <- read.table(samples_file, head=FALSE, stringsAsFactors=FALSE)
metadata_sampled <- metadata[metadata$Accession.ID %in% samples$V1,]
#stopifnot(dim(samples)[1] == dim(metadata_sampled)[1])
write.table(metadata_sampled, metadata_sampled_file, sep="\t", quote=FALSE, row.names = FALSE)

