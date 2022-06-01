

#' Plot Cap-TCRseq QC metrics
#' @description This function generates series of plots as QC metrics for Cap-TCRseq runs. Input files are csv files generated from log_parser.py.
#' Samples are ordered by total sequencing reads.
#'
#' 
#' @param datapath path to log_parser output files
#' @param alignstatsfile align_stats.csv filename (from log_parser)
#' @param assemblestatsfile assemble_stats.csv filename (from log_parser)
#' @param plotname name for output plot
#' @param plotpath path to plot directory
#'
#'
#' @examples

mixcrQC.fx <- function(datapath, alignstatsfile, assemblestatsfile, plotname, plotpath){
  
  if (!file.exists(paste0(plotpath, plotname))){
    
    # Function to align plots (from stackoverflow) 
    align_plots1 <- function (...) {
      pl <- list(...)
      stopifnot(do.call(all, lapply(pl, inherits, "gg")))
      gl <- lapply(pl, ggplotGrob)
      bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
      combined <- Reduce(bind2, gl[-1], gl[[1]])
      wl <- lapply(gl, "[[", "widths")
      combined$widths <- do.call(grid::unit.pmax, wl)
      grid::grid.newpage()
      grid::grid.draw(combined)
    }
    
    mytheme <- theme(axis.title.y = element_text(size = 25),
                     axis.title.x = element_blank(),
                     axis.line = element_line(color = "black"),
                     axis.text = element_text(size = 22),
                     axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),  
            legend.key = element_rect(fill = "white", colour = "white"),
            plot.margin = unit(c(0,0,0,0),"cm"))    
    
    #read stats files
    alignstats <- read.csv(file = paste0(datapath,alignstatsfile),
                           sep = ",", header = T, stringsAsFactors = F)    
    assemblestats <- read.csv(file = paste0(datapath, assemblestatsfile),
                              sep = ",", header = T, stringsAsFactors = F)
    
    #Create samplenames
    alignstats$samplename <- gsub(".*alignments_", "", alignstats$Output.file.s.)
    alignstats$samplename <- gsub(".vdjca", "", alignstats$samplename) 
    
    assemblestats$samplename <- gsub(".*extended_", "", assemblestats$Input.file.s.)
    assemblestats$samplename <- gsub(".vdjca", "", assemblestats$samplename) 
    
    #Order samplenames based on total seq reads
    alignstats$samplename <- factor(alignstats$samplename, 
                                    levels = alignstats$samplename[order(alignstats$Total.sequencing.reads)])
    assemblestats$samplename <- factor(assemblestats$samplename, 
                                       levels = alignstats$samplename[order(alignstats$Total.sequencing.reads)])           
    
    #plot total sequencing reads plot
    myplot_totalseq <- ggplot(data = alignstats, aes(x = samplename)) +
      geom_point(aes(y = Total.sequencing.reads), color = "black",size = 5) + 
      mytheme + theme(axis.text.x = element_blank()) +
      # This is to avoid printing min, max and median on y-axis
      scale_y_continuous(limits = c(0, max(alignstats$Total.sequencing.reads)),
                         sec.axis = dup_axis(name = "Aligned reads"))    
    
    
    # barplot percentage aligned reads
    myplot_percaligned <- ggplot(aes(x = samplename, 
                                     y = Successfully.aligned.reads/Total.sequencing.reads), 
                                 data = alignstats) + 
      geom_bar(stat = "identity") + 
      mytheme + theme(axis.text.x = element_blank())        
    
    
    #plot average reads per clonotype  
    assemblestats$Average.number.of.reads.per.clonotype <- as.numeric(assemblestats$Average.number.of.reads.per.clonotype)  
    myplot_averagereads <- ggplot(aes(x = samplename, y = Average.number.of.reads.per.clonotype), 
                                  data = assemblestats) + 
      geom_point(size = 5) + 
      mytheme + theme(axis.text.x = element_blank())
    
    #plot reads used in clonotypes
    myplot_readsused <- ggplot(aes(x = samplename, y = Reads.used.in.clonotypes.before.clustering..percent.of.total),       
                               data = assemblestats) + 
      geom_point(size = 5) + 
      mytheme + theme(axis.text.x = element_blank())
    
    #plot total clonotypes
    myplot_totalclon <- ggplot(aes(x = samplename, y = Final.clonotype.count),        
                               data = assemblestats) + 
      geom_point(size = 5) + 
      mytheme
    
    
    pdf(file = paste(plotpath, plotname, sep = ""),
        width = 15, 
        height = 20,
        useDingbats = FALSE,
        onefile = FALSE)
    align_plots1(myplot_totalseq + #Add aligned reads as secondary y axis
                   geom_point(aes(y = Successfully.aligned.reads), color = "red",size = 5) + 
                   theme(axis.title.y.right = element_text(color = "red"))+
                   ylab("Total sequencing reads"),
                 myplot_percaligned + ylab("Percentage\naligned reads"),
                 myplot_readsused + ylab("Reads used \nin clonotypes"),
                 myplot_averagereads + ylab("Average \nreads/clonotype"),
                 myplot_totalclon + ylab("Total clonotypes"))    
    dev.off()     
  }  else {
    message("plot file exists: ", plotname)
  }
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
logdir =args[1]
mixcrQC.fx(logdir, "align_stats.csv", 
           "assemble_stats.csv", 
           "project_QC.pdf", logdir)


