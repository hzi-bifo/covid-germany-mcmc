require(rjson)
require(treeio)
require(ggplot2)
require(ggsci)
require(ggtree)
require(phytools)

load_metadata <- function(filename) {
  ############
  # Metadata #
  ############
  #metadata              <- read.csv(paste0(inputpath, "metadata.csv"))
  metadata <- read.table(filename, sep="\t", head=TRUE, na.strings=c("NA", ""), fill=TRUE, stringsAsFactors=FALSE, quote="|")
  metadata <- metadata[metadata$Virus.name != "Virus.name",]
  metadata$sample_date_orig  <- ymd(metadata$Collection.date)
  # metadata$sample_date  <- ymd(metadata$date_corrected)
  metadata$sample_date  <- metadata$sample_date_orig
  metadata$decimal_date <- decimal_date(metadata$sample_date)    
  metadata$taxon_label  <- metadata$Accession.ID
  metadata$country  <- str_trim(sapply(str_split(metadata$Location, "/"), "[[", 2))
  # metadata$state    <- str_trim(sapply(str_split(paste0(metadata$Location,"/",metadata$Additional.location.information), "/"), "[[", 3))
  metadata$adm1 <- sapply(
    str_split(paste0(metadata$Location,"/",metadata$Additional.location.information), "/"), function(x) str_trim(x[3]) )
  
  return(metadata);
}

set_instate <- function(metadata, state) {
  return(sapply(str_split(paste0(metadata$Location,"/",metadata$Additional.location.information), "/"), function(x) {return(str_trim(x[2]) == "Germany" & (grepl(state, x[3], fixed=TRUE) | ( length(x) >= 4 & grepl(state, x[4], fixed=TRUE) )) );}))
}
# metadata$instate      <- set_instate(metadata, state)

load_case_data <- function(filename) {
  case_data <- read.table(filename, sep = ",", quote = "\"", stringsAsFactors = FALSE)
  case_data_date <- as.Date(unlist(case_data[1,5:ncol(case_data)]), format="%m/%d/%y")
  case_data_value <- sapply(case_data[2:nrow(case_data), 5:ncol(case_data)], as.numeric)
  case_data_country <- case_data[2:nrow(case_data),2]
  case_data <- aggregate(case_data_value, data.frame(country = case_data_country), sum)
  case_data_country <- case_data[,1]
  case_data <- case_data[,2:ncol(case_data)]
  case_data[,2:ncol(case_data)] <- case_data[,2:ncol(case_data)] - case_data[,1:(ncol(case_data)-1)]
  rownames(case_data) <- case_data_country
  return(list(data=case_data, country=case_data_country, date = case_data_date));
}

average_ <- function(v_, len = 7) {
  v <- c(v_, rep(0, len-1))
  v_s <- cumsum(v)
  v_s_bet <- v_s
  v_s_bet[(len+1):length(v)] <- v_s[(len+1):length(v)] - v_s[1:(length(v) - len)]

  l <- c(rep(1, length(v_)), rep(0, len-1))
  l_s <- cumsum(l)
  l_s_bet <- l_s
  l_s_bet[(len+1):length(l)] <- l_s[(len+1):length(l)] - l_s[1:(length(l) - len)]
  
  r <- v_s_bet / l_s_bet
  return(r[((len+1)%/%2):(length(v_) + (len+1)%/%2-1)])
}

plot_on_map <- function(names, values, title, legend_title, col = "#101010", min_value = 0, max_value = NULL, legend.position = "right", text.font.size = NULL) {
  germany <- raster::getData("GADM", country = "DEU", level = 1)
  germany.f <- fortify(germany, region = "CC_1")
  
  names[names == "Baden-Württemberg"] <- "Baden-Wurttemberg"
  
  mapNames <- germany$VARNAME_1
  mapNames[is.na(mapNames)] <- germany$NAME_1[is.na(mapNames)]
  mapNames[mapNames == "Baden-Württemberg"] <- "Baden-Wurttemberg"
  mapNames[mapNames == "Mecklenburg-West Pomerania"] <- "Mecklenburg-Western Pomerania"
  
  germany$NAME_1_EN <- germany$NAME_1
  germany$NAME_1_EN[germany$NAME_1 == "Baden-Württemberg"] <- "Baden-Wurttemberg"
  germany$NAME_1_EN[germany$NAME_1 == "Mecklenburg-Vorpommern"] <- "Mecklenburg-Western Pomerania"
  germany$NAME_1_EN[germany$NAME_1 == "Bayern"] <- "Bavaria"
  germany$NAME_1_EN[germany$NAME_1 == "Hessen"] <- "Hesse"
  germany$NAME_1_EN[germany$NAME_1 == "Niedersachsen"] <- "Lower Saxony"
  germany$NAME_1_EN[germany$NAME_1 == "Nordrhein-Westfalen"] <- "North Rhine-Westphalia"
  germany$NAME_1_EN[germany$NAME_1 == "Rheinland-Pfalz"] <- "Rhineland-Palatinate"
  germany$NAME_1_EN[germany$NAME_1 == "Sachsen"] <- "Saxony"
  germany$NAME_1_EN[germany$NAME_1 == "Sachsen-Anhalt"] <- "Saxony-Anhalt"
  germany$NAME_1_EN[germany$NAME_1 == "Thüringen"] <- "Thuringia"
  
  germany$value <- 0
  germany$value[match(names, germany$NAME_1_EN)] <- values
  
  
  germany.f$value <- germany$value[match(germany.f$id, germany$CC_1)]
  if (is.null(max_value)) {
    max_value <- max(germany$value) * 1.1
  }
  
  # id=11 is drawn again because it is inside another region.
  p <- ggplot() + 
    geom_polygon(data = germany.f, aes(x = long, y = lat, group = group, fill = value), colour = "grey10")+
    geom_polygon(data = germany.f[germany.f$id == "11",], aes(x = long, y = lat, group = group, fill = value), colour = "grey10")+
    scale_fill_gradient(low = lighten(col, 1), high = col, space = "Lab", limits=c(min_value, max_value),
                        name = legend_title)+
    ggtitle(title) +
    # options(repr.plot.width=max(germany.f$lat) - min(germany.f$lat), repr.plot.height=max(germany.f$long) - min(germany.f$long)) +
    theme(axis.line=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(), panel.grid.major=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.1),
          panel.grid.minor=element_blank(), plot.background=element_blank(), 
          aspect.ratio = (max(germany.f$long) - min(germany.f$long))/(max(germany.f$lat) - min(germany.f$lat)),
          legend.position = legend.position,
          text = element_text(size = text.font.size)
          )
  print(p)
  
}

getDeInfo <- function(populationFileName, deCountyDaily, startDate) {
  dePopulation <- read.csv(populationFileName)
  
  deCounty <- deCountyDaily %>% filter(date >= startDate) %>% group_by(county) %>% dplyr::summarise(case = sum(case), seqs = sum(seqs), seq.rate = sum(seqs) / sum(case)) %>% dplyr::rename(state = county) %>% filter(!is.na(state))
  
  
  deCounty <- deCounty %>% left_join(dePopulation %>% dplyr::rename(population = X2020), by=c("state" = "State"))
  
  return(deCounty)
}

getDeCountyDaily <- function(caseDataGermany, germanySubdivisions, metadata) {
  case_data_germany <- read.table(caseDataGermany, sep=",", quote = "\"", header = TRUE)
  case_data_germany[2:nrow(case_data_germany),2:ncol(case_data_germany)] <- case_data_germany[2:nrow(case_data_germany),2:ncol(case_data_germany)] - case_data_germany[1:(nrow(case_data_germany)-1),2:ncol(case_data_germany)]
  
  deCountyDaily <- case_data_germany %>% gather(key=county, value=case, -time_iso8601) %>% mutate(date = as.Date(time_iso8601)) %>% dplyr::select(-time_iso8601)
  
  deCountyNames <- read.csv(germanySubdivisions)
  deCountyDaily$county <- deCountyNames$subdivision.name.en[match(str_replace(deCountyDaily$county, "\\.", "-"), deCountyNames$adm1)]
  
  deCountyDaily <- deCountyDaily %>% left_join(metadata %>% filter(country == "Germany") %>% group_by(adm1, sample_date) %>% dplyr::summarise(seqs = n()), by=c("county" = "adm1", "date" = "sample_date")) %>% mutate(seqs = ifelse(is.na(seqs), 0, seqs))
  
  return(deCountyDaily);
}

# loadMetadata <- function(filename) {
#   #metadata              <- read.csv(paste0(inputpath, "metadata.csv"))
#   metadata <- read.table(filename, sep="\t", head=TRUE, na.strings=c("NA", ""), fill=TRUE, stringsAsFactors=FALSE, quote="|")
#   metadata <- metadata[metadata$Virus.name != "Virus.name",]
#   metadata$sample_date  <- ymd(metadata$Collection.date)
#   metadata$decimal_date <- decimal_date(metadata$sample_date)    
#   #metadata$taxon_label  <- metadata$sequence_name
#   metadata$taxon_label  <- metadata$Accession.ID
#   #metadata$taxon_label <- gsub("/", "_", as.character(metadata$sequence_name))
#   metadata$country                 <- sapply(strsplit(metadata$Location, '/'), function(x) str_trim(x[2]))
#   # metadata$state      <- str_trim(sapply(str_split(paste0(metadata$Location,"/",metadata$Additional.location.information), "/"), "[[", 3))
#   return(metadata) 
# }


calculateWeek <- function(dates, referenceDate) {
  return(as.numeric(dates - referenceDate)%/%7)
}

