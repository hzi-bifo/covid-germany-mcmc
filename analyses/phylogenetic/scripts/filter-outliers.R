suppressWarnings(suppressMessages(library(dplyr)))
args = commandArgs(trailingOnly=TRUE)
lineage_info_file = args[1]
metadata_file = args[2]
metadata_out_file = args[3]

lineage_info <- read.table(lineage_info_file, sep='|', quote='"', fill=TRUE, head=TRUE, stringsAsFactors=FALSE)
metadata <- read.table(metadata_file, sep="\t", head=TRUE, na.strings=c("NA", ""), fill=TRUE, stringsAsFactors=FALSE, quote="|")

#earliest <- min(lineage_info[startsWith(lineage_info$Lineage, paste0(lineage_name,".")) | lineage_info$Lineage == lineage_name,]$Earliest.date)
#if (sum(metadata$Collection.date < earliest) > 0) {cat(paste("Removed", sum(metadata$Collection.date < earliest), "samples from", lineage_name, "\n"), file=stderr()); }
#write.table(metadata[metadata$Collection.date>=earliest,], metadata_out_file, sep="\t", quote=FALSE, row.names = FALSE)

earliest <- sapply(levels(as.factor(metadata$Pango.lineage)), function(lineage_name) {suppressWarnings(earliest <- min(lineage_info[startsWith(lineage_info$Lineage, paste0(lineage_name,".")) | lineage_info$Lineage == lineage_name,]$Earliest.date)); return(earliest);} )
earliest2 <- sapply(names(earliest[is.na(earliest)]), function(lineage_name) { while (lineage_name != "" && ( !(lineage_name %in% names(earliest)) || (is.na(earliest[lineage_name])) ) ) { if (grepl("\\.", lineage_name)) { lineage_name = substr(lineage_name, 1, regexpr("\\.[^.]*$", lineage_name)[1]-1);} else {lineage_name = "";} }; if (lineage_name %in% names(earliest)) return(earliest[lineage_name][[1]]) else return(NA); } )
for (n in names(earliest2)) {earliest[n] = earliest2[n];}
d <- sapply(names(earliest[earliest == ""])[!is.na(names(earliest[earliest == ""]))], function(x) ifelse(lineage_info[lineage_info$Lineage == x, 'Earliest.date'] != "", lineage_info[lineage_info$Lineage == x, 'Earliest.date'], ""))
earliest[names(d)] = d

earliest = unlist(earliest)

earliest.df = data.frame(lineage=names(earliest), earliest=earliest, stringsAsFactors = FALSE)

metadata_filtered <- metadata %>% left_join(earliest.df, by=c('Pango.lineage' = 'lineage')) %>% filter(earliest == "" | is.na(earliest) | Collection.date >= earliest) %>% select(-earliest)
if (nrow(metadata) != nrow(metadata_filtered)) {
	cat(paste("Removed", nrow(metadata) - nrow(metadata_filtered), "samples from", metadata_file, "\n"), file=stderr());
}
write.table(metadata_filtered, metadata_out_file, sep="\t", quote=FALSE, row.names = FALSE)

