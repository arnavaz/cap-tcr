library(ggplot2)
require(ggalluvial)
library(randomcoloR)


#' Dataframe for clonetracking
#' @description This function takes a list of files and compiles them in one dataframe
#' for use in plot_clonetracks.fx. Also you can choose to exclude non-productive clonotypes
#'
#' @param datapath path to mixcr files
#' @param chain any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#' @param filelist a list of files to track clonotypes. Output from list.files()
#' @param totalinframe "total" outputs all CDR3s, "inframe" removes non-productive CDR3s
#'
#'
#' @examples

cdr3_dataframe.fx <- function(datapath, chain, filelist, totalinframe){

  if (!(totalinframe %in% c("total", "inframe"))) {
    stop("Error: unknown argument ", totalinframe, ". Please provide either total (for all clonotypes) or inframe (for in-frame clonotypes only)")
  }

  # Ensure only one chain is included
  filelist <- filelist[grepl(chain, filelist)]

  #Compile a big file with patient's mixcr files
  i <- 1
  for (f in filelist){
    mixcrfle <- read.table(paste(datapath, f, sep = ""),
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))
    if(i == 1){
      compldfle <- mixcrfle
      compldfle <- cbind(cloneno = row.names(compldfle), filename = f, compldfle)
      i <- i + 1
    }
    else{
      compldfle1 <- mixcrfle
      compldfle1 <- cbind(cloneno = row.names(compldfle1), filename = f, compldfle1)
      compldfle <- rbind(compldfle, compldfle1)
      rm(compldfle1)
    }
  }
  myfiles <- unique(as.character(compldfle$filename))
  message("my files:")
  print(myfiles)

  message("Total recovered clonotypes:")
  print(length(compldfle$aaSeqCDR3))

  message("Total out-of-frame clonotypes:")
  print(length(compldfle$aaSeqCDR3[grepl("_", compldfle$aaSeqCDR3)]))
  message("Total clonotypes with stop codon:")
  print(length(compldfle$aaSeqCDR3[grepl("[*]", compldfle$aaSeqCDR3) &
                                     !grepl("_", compldfle$aaSeqCDR3)]))

  #make samplename column
  compldfle$filename <- as.character(compldfle$filename)
  compldfle$samplename <- gsub(".*.CLONES_","", compldfle$filename)

  # remove out-of-frame clonotypes and those with stop codon
  compldfle_clean <- compldfle[!grepl("_", compldfle$aaSeqCDR3) &
                                 !grepl("[*]", compldfle$aaSeqCDR3),]
  #Recalculate cloneFraction for each file
  compldfle_clean$cloneFraction <- NA
  for(f in myfiles){
    compldfle_clean$cloneFraction[compldfle_clean$filename == f] <- compldfle_clean$cloneCount[compldfle_clean$filename == f]/sum(compldfle_clean$cloneCount[compldfle_clean$filename == f])
  }

  message("Total productive clonotypes:")
  print(length(compldfle_clean$aaSeqCDR3))

  if(totalinframe == "inframe"){
    message("Output contains in_frame clonotypes only")
    return(compldfle_clean)}
  if(totalinframe == "total"){
    message("Output contains all clonotypes")
    return(compldfle)}
}


#' Plot recurring CDR3s
#' @description Hello!
#'
#'
#' This function uses the dataframe from cdr3_dataframe.fx function,
#' identifies recurring nucleotideCDR3s (found in at least two samples) and color them on clones stacked barplot. If no recurring CDR3s are found,
#' it outputs the clones stacked barplot with no colors. If recurring CDR3 are more than 50, it tracks top 10 clonotypes. It recalculates cloneFraction based on cloneCount column in MiXCR output.
#'
#' If you need to clean up samplenames or change the order of samples on the stacked barplot,
#' just modify your dataframe.
#'
#' @param compldfle dataframe from cdr3_dataframe.fx
#' @param plotpath path to plot directory
#' @param chain any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#' @param countfrac plot either clonal fraction of absolute counts
#' ("cloneCount" or "cloneFraction")
#' @param clnefrc specify a cut-off from 0 to 1 to track and plot only a subset of clonotypes.
#' Useful when you have too many clonotypes to plot.
#' Clonal fraction of 0.001 is usually a good starting point. If you want all just pass 0.
#'
#'
#' @examples

plot_clonetracks.fx <- function(compldfle, plotpath, chain, countfrac, clnefrc){

  if (!(countfrac %in% c("cloneFraction", "cloneCount"))) {
    stop("Error: unknown argument ", countfrac, ". Please provide either cloneFraction or cloneCount.")
  }

  message("list of samples to track clones: ")
  mysamples <- unique(compldfle$samplename)
  print(mysamples)

  # Subset df
  CDR3_fraction <- compldfle[, c("samplename","nSeqCDR3","cloneFraction", "cloneCount")]

  # Subset to include only clonotypes with more than specified clonal fraction
  CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > clnefrc,]

  #Assign colors to recurring nt clonotypes
  recurring <- unique(CDR3_fraction$nSeqCDR3[duplicated(CDR3_fraction$nSeqCDR3)])
  notrecurring <- CDR3_fraction$nSeqCDR3[!CDR3_fraction$nSeqCDR3 %in% recurring]

  message("Total number of recurring clonotypes: ")
  print(length(recurring))

  if(length(recurring) == 0){
    #Introduce a dummy common cdr3 dataframe for alluvia
    mydummy_df <- as.data.frame(matrix(ncol = 4, nrow = length(mysamples)))
    colnames(mydummy_df) <- colnames(CDR3_fraction)

    mydummy_df$samplename <- mysamples
    mydummy_df$nSeqCDR3 <-  "XXXXX"
    mydummy_df$cloneFraction <- 0
    mydummy_df$cloneCount <- 0
    CDR3_fraction <- rbind(CDR3_fraction, mydummy_df)

    recurring <- "XXXXX"
  }

  if(length(recurring) > 50){
    recurring_df <- CDR3_fraction[CDR3_fraction$nSeqCDR3 %in% recurring,]
    recurringcdr3_ordered <- unique(recurring_df$nSeqCDR3[order(recurring_df$cloneCount, decreasing = TRUE)])
    message("Total number of recurring clonotypes > 50 ")
    message("Tracking top 10 recurring clonotypes ")
    myColors <- distinctColorPalette(10)

    myColors <- c(myColors, rep("white",length(recurring)-10),
                  rep("white",length(notrecurring)))
    names(myColors) <- c(recurringcdr3_ordered, notrecurring)

    message("these are what we color: ")
    print(myColors[myColors != "white"])
  }
  else{
    myColors <- distinctColorPalette(length(recurring))
    myColors <- c(myColors, rep("white",length(notrecurring)))
    names(myColors) <- c(recurring, notrecurring)

    myColors[names(myColors) == "XXXXX"] <- "white"

    message("these are what we color: ")
    print(myColors[myColors != "white"])
  }

  # Generate a row for each sample that doesnot have recurring clonotype
  ## This ensures alluvia are colored

  for(c in recurring){
    tmp <- CDR3_fraction[CDR3_fraction$nSeqCDR3 == c,]
    nonexsiting <- mysamples[!mysamples %in% tmp$samplename]
    if(length(nonexsiting) > 0){
      newentries <- data.frame("samplename" = nonexsiting, "nSeqCDR3" = c,
                               "cloneFraction" = 0, "cloneCount" = 0)
      CDR3_fraction <- rbind(CDR3_fraction, newentries)
    }
  }


  p <-  ggplot(CDR3_fraction, aes(x = samplename,
                                  y = eval(as.name(countfrac)),
                                  fill = nSeqCDR3,
                                  stratum = nSeqCDR3,
                                  alluvium = nSeqCDR3,
                                  label = nSeqCDR3))

  myp <- p + geom_alluvium(decreasing = FALSE) +
    geom_stratum(decreasing = FALSE, stat = "alluvium") +
    scale_fill_manual(breaks = names(myColors[myColors != "white"]),
                      values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "none",
          plot.margin = unit(c(0.2,0,0,0),"cm")) +
    labs(y = countfrac)

  pdf(paste0(plotpath, "clonetrack_", mysamples[1],
             chain, countfrac, ".pdf"),
      width = 15,
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)
  print(myp)
  dev.off()

}

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
datapath=args[1]
pattern=args[2]
type=args[3]
flelst <- list.files(datapath,
                     recursive = F,
                     pattern = paste0(type, ""))

pt_df <- cdr3_dataframe.fx(datapath, type, flelst[grepl(patt, flelst)], "inframe")
plot_clonetracks.fx(pt_df, datapath, type, "cloneFraction", 0) 

