################################################################################################################################
# Palettes
################################################################################################################################

dark  <- list(blue   = RColorBrewer::brewer.pal(12,"Paired")[2], 
              green  = RColorBrewer::brewer.pal(12,"Paired")[4], 
              red    = RColorBrewer::brewer.pal(12,"Paired")[6], 
              orange = RColorBrewer::brewer.pal(12,"Paired")[8], 
              purple = RColorBrewer::brewer.pal(12,"Paired")[10], 
              gray   = "#777777",
              black  = "#000000",
              white  = "#FFFFFF")

light <- list(blue   = RColorBrewer::brewer.pal(12,"Paired")[1], 
              green  = RColorBrewer::brewer.pal(12,"Paired")[3], 
              red    = RColorBrewer::brewer.pal(12,"Paired")[5], 
              orange = RColorBrewer::brewer.pal(12,"Paired")[7], 
              purple = RColorBrewer::brewer.pal(12,"Paired")[9], 
              gray   = "#777777",
              black  = "#000000",
              white  = "#FFFFFF")


ukPal <- list(eng = "#BE0F34", 
              sct = "#191970", 
              wls = "#F5CF47",
              nir = "#9ECEEB", 
              oth = "#C7C2BC")

countryPal <- list("China"         = "#d53e4f",
                   "Italy"         = "#fc8d59",
                   "Spain"         = "#fee08b",  
                   "France"        = "#99d594",
                   "Belgium"       = "#3288bd",
                   "Netherlands"   = "#FF7F00",
                   "Ireland"       = "#AAB300",
                   "Switzerland"   = "#BE0F34",
                   "Germany"       = "#CF7A30", 
                   "US"            = "#A6CEE3",
                   "Sweden"        = "#007770",
                   "Portugal"      = "#6A3D9A",
                   "Other"         = "#C7C2BC")

dePal <- list(de = "#fde0ef",
              "Germany" = "#fde0ef",
              
              "North_Rhine-Westphalia" = "#c51b7d",
              "North Rhine-Westphalia" = "#c51b7d",    
              "Baden-Wurttemberg" = "#e9a3c9",
              "Baden-Württemberg" = "#e9a3c9",             
             
              "Saarland" = "#8e0152",
              "Rhineland-Palatinate" = "#de77ae",
              "Lower_Saxony" = "#bda6b2",
              "Lower Saxony" = "#bda6b2",
              "Hesse" = "#fa9fb5",
              "Bremen" = "#7a0177",
              "Hamburg" = "#7a0177",
              
              "Bavaria" = "#4d9221",
              "Saxony" = "#00441b",
              "Brandenburg" = "#b8e186",
              "Mecklenburg-Vorpommern" = "#cbe6a5",
              "Mecklenburg-Western Pomerania" = "#cbe6a5",
              "Saxony-Anhalt" = "#a1d76a",
              "Schleswig-Holstein" = "#737d70",
              "Thuringia" = "#7fbc41",
              "Berlin" = "#276419",
              
              oth = "#C0C0C0",
              all = "#1d1d1d",
              "NA" = "#252525",
              "Non Bavaria" = "#C0C0C0",
              "Dusseldorf"  = "#FFb266",
              "Non Dusseldorf" = "#C0C0C0",
              
              "Non Germany" = "#B0B0B0",
              
              "Non Hamburg" = "#C0C0C0",
              
              "Non Lower_Saxony" = "#C0C0C0",
              
              "Non Lower Saxony" = "#C0C0C0",
              "Munich"  = "#66FF66",
              "Non Munich" = "#C0C0C0",
              "Non North_Rhine-Westphalia" = "#C0C0C0",
              "Non North Rhine-Westphalia" = "#C0C0C0",
              "Non Saarland" = "#C0C0C0",
              "Non Baden-Württemberg" = "#C0C0C0",
              "Non Baden-Wurttemberg" = "#C0C0C0",
              "Non Berlin" = "#C0C0C0",
              "Non Brandenburg" = "#C0C0C0",
              "Non Bremen" = "#C0C0C0",
              "Non Hesse" = "#C0C0C0",
              "Non Mecklenburg-Vorpommern" = "#C0C0C0",
              "Non Mecklenburg-Western Pomerania" = "#C0C0C0",
              "Non Rhineland-Palatinate" = "#C0C0C0",
              "Non Saxony" = "#C0C0C0",
              "Non Saxony-Anhalt" = "#C0C0C0",
              "Non Schleswig-Holstein" = "#C0C0C0",
              "Non Thuringia" = "#C0C0C0",
              "Herzogenrath" = "#A0C0B0"
)

lineagePal <- list(
              "B.1.1.7_MCC_237" = "#d53e4f",
              "B.1.1.7_MCC_164" = "#fc8d59",    
              "B.1.1.7_MCC_181" = "#fee08b",
              "B.1.1.7_MCC_116" = "#99d594",             
              "B.1.1.7_MCC_16" = "#3288bd",
              "B.1.1.7_MCC_193" =  "#5e4fa2",
              "A_MCC_2" = "#3288bd",
              "Other" = "#C0C0C0",
              "NA" = "#C0C0C0",
              "B.1.1.7_MCC_363" = "#ffffff",
              "B.1.1.7_MCC_2" = "#f0f0f0",
              "B.1.1.7_MCC_384" = "#d9d9d9",
              "B.1.1.7_MCC_48" = "#bdbdbd",
              "B.1.1.7_MCC_86" = "#969696",
              "B.1.1.7_MCC_143" = "#737373",
              "B.1.1.7_MCC_301" = "#525252"
             
              
             
              
)


pangoPal <- list(
  "C"         = "#fcd933",
  "A"         = "#fffcd6",
  "B.1.1.70"  = "#ed93C9",
  
  "B.1.1.317" = "#80b1d3",
  "B.1.36"    = "#fdb462",  
  
  "B.1.177"   = "#d44076",
  "B.1.221"   = "#f6ccd7",
  "B.1.1.7"   = "#eb4ba0",
  
  "B.1.160"   = "#37a1fb",
  "B.1.258"   = "#154e80",
  "B.1.1.519" = "#d9d9d9",
  "B.1.351"   = "#2676ba",
  "Other"     = "#bdbdbd"
)

sdePal <- list(de = "#c51b7d",
              "Germany" = "#c51b7d",
              
              "North_Rhine-Westphalia" = "#e9a3c9",
              "North Rhine-Westphalia" = "#e9a3c9",    
              "Baden-Wurttemberg" = "#a1d76a",
              "Baden-Württemberg" = "#a1d76a",             
              
              "Saarland" = "#bdbdbd",
              "Rhineland-Palatinate" = "#bdbdbd",
              "Lower_Saxony" = "#bdbdbd",
              "Lower Saxony" = "#bdbdbd",
              "Hesse" = "#bdbdbd",
              "Bremen" = "#bdbdbd",
              "Hamburg" = "#bdbdbd",
              
              "Bavaria" = "#4d9221",
              "Saxony" = "#bdbdbd",
              "Brandenburg" = "#bdbdbd",
              "Mecklenburg-Vorpommern" = "#bdbdbd",
              "Mecklenburg-Western Pomerania" = "#bdbdbd",
              "Saxony-Anhalt" = "#bdbdbd",
              "Schleswig-Holstein" = "#bdbdbd",
              "Thuringia" = "#bdbdbd",
              "Berlin" = "#bdbdbd",
              
              oth = "#bdbdbd",
              "Other" = "#bdbdbd",
              all = "#bdbdbd",
              "NA" = "#bdbdbd",
              "Non Bavaria" = "#bdbdbd",
              "Dusseldorf"  = "#bdbdbd",
              "Non Dusseldorf" = "#bdbdbd",
              
              "Non Germany" = "#bdbdbd",
              
              "Non Hamburg" = "#bdbdbd",
              
              "Non Lower_Saxony" = "#bdbdbd",
              
              "Non Lower Saxony" = "#bdbdbd",
              "Munich"  = "#bdbdbd",
              "Non Munich" = "#bdbdbd",
              "Non North_Rhine-Westphalia" = "#bdbdbd",
              "Non North Rhine-Westphalia" = "#bdbdbd",
              "Non Saarland" = "#bdbdbd",
              "Non Baden-Württemberg" = "#bdbdbd",
              "Non Baden-Wurttemberg" = "#bdbdbd",
              "Non Berlin" = "#bdbdbd",
              "Non Brandenburg" = "#bdbdbd",
              "Non Bremen" = "#bdbdbd",
              "Non Hesse" = "#bdbdbd",
              "Non Mecklenburg-Vorpommern" = "#bdbdbd",
              "Non Mecklenburg-Western Pomerania" = "#bdbdbd",
              "Non Rhineland-Palatinate" = "#bdbdbd",
              "Non Saxony" = "#bdbdbd",
              "Non Saxony-Anhalt" = "#bdbdbd",
              "Non Schleswig-Holstein" = "#bdbdbd",
              "Non Thuringia" = "#bdbdbd",
              "Herzogenrath" = "#bdbdbd"
)
################################################################################################################################

mPal <- function(c, alpha=1.0) {
  if (length(c) > 1) {
    return(sapply(c, function(x) mPal(x, alpha)))
  } else if (is.character(c) && substr(c,1,1) == "#") {
      return(paste0(c,format(as.hexmode(round(alpha*255)), width=2)))
  } else {
      return(rgb(red=c[1], green=c[2], blue=c[3], alpha=round(alpha*255), maxColorValue=255))
  }
}


plotPalette <- function(pal, alpha=1.0) {
  
  root <- sqrt(length(pal))
  layout(matrix(1:(round(root)*ceiling(root)), nrow=round(root)))
  
  par(mar=c(2,0,0,0))
  for (col in 1:length(pal)) {
      plot(1,type='n',xlim=c(0,1),ylim=c(0,1), axes=FALSE, ylab="", xlab="")
      rect(0,0,1,1,col=mPal(pal[[col]], alpha=alpha))
      if (is.null(names(pal)[col])) {
          mtext(col,line=0,side=1)
      } else {
          mtext(names(pal)[col],line=0,side=1)
      }
  }
}


