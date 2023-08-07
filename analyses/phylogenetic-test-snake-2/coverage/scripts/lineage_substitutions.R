#file <- "/home/sreimering/Documents/Coronavirus/Hadi/importation_lineages_DE_casesampling/data/clusterSamples_DTA_MCC_0.5_reformatted.tsv"
file <- "data/clusterSamples_DTA_MCC_0.5_reformatted.tsv"

# read data
gisaid_data <- read.csv2(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# function to get substitutions within a lineage
substitution_count <- function(gisaid_data, substitution){
  substitutions <- gisaid_data[gisaid_data$pangolin_lineage == substitution,]$AA.Substitutions
  substitutions <- gsub("(", "", substitutions, fixed = TRUE)
  substitutions <- gsub(")", "", substitutions, fixed = TRUE)
  
  substitutions_list <- strsplit(substitutions, split = ",", fixed = TRUE)
  all_substitutions <- unlist(substitutions_list)
  substitution_counts <- table(all_substitutions)[order(table(all_substitutions))]
  
  # substitutions shared in a lineage occur in > 90% of sequences
  shared_subs <- names(substitution_counts)[substitution_counts/length(substitutions) > 0.9]
  
  # investigate if these substitutions also occur in other lineages
  substitutions_other <- gisaid_data[gisaid_data$pangolin_lineage != substitution,]$AA.Substitutions
  substitutions_other <- gsub("(", "", substitutions_other, fixed = TRUE)
  substitutions_other <- gsub(")", "", substitutions_other, fixed = TRUE)
  substitutions_list_other <- strsplit(substitutions_other, split = ",", fixed = TRUE)
  all_substitutions_other <- unlist(substitutions_list_other)
  
  count_other <- vector(mode = "integer", length=length(shared_subs))
  for (i in 1:length(shared_subs)){
    count_other[i] <- sum(all_substitutions_other == shared_subs[i])
  }
  
  results <- cbind(shared_subs, count = unname(substitution_counts[shared_subs]), frac = unname(substitution_counts[shared_subs]/(substitution_counts[shared_subs] + count_other)), count_other)
  return(results)
}


### significant lineages in the new analysis

lineages <- c("B.1.1.7_237", "B.1.1.7_363", "B.1.1.7_116", "B.1.1.7_164")

# this outputs onto the console:
# 1. the name of the importation lineage
# 2. a table with all shared substitutions, including the number of sequences having this substitutions within the lineage (count), 
#    the number of sequences having this substitution outside the lineage (count_other), and
#    the fraction of sequences with that substituion that are in the lineage (fraction 1 means this substitutions is lineage-specific)
# 3. a string containing all spike substitutions
for (lin in lineages){
  table(gisaid_data[gisaid_data$pangolin_lineage == lin,]$Pango.lineage)
  results <- substitution_count(gisaid_data, lin)
  print(lin)
  print(results)
  print(paste(gsub("Spike_", "", results[startsWith(results[,1], "Spike_"),1], fixed = TRUE), collapse = ", "))
}
