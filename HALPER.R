#!/usr/bin/env Rscript
# implementes HALPER method for extending hal mapped summits

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))

option_list <- list(
  make_option(c("-q","--qFile"), 
              help="input bed file, 1st 4 columns must be in standard bed format, peak names must be unique"),
  make_option(c("-t","--tFile"), help = "input mapped bed file"),
  make_option(c("-s","--sFile"), help = "input mapped-summit bed file"),
  make_option(c("-o","--oFile"), help = "output file"),
  make_option("--max_len",default = 1500, help = " maximum number of base pairs of the ortholog. if < 10, treated as max fraction. Default 1500"),
  make_option("--min_len", default = 100, help = "minimum number of base pairs of the ortholog; if < 1, treated as min fraction. Default 100"),
  make_option("--protect_dist", default = 25, help = "summit protection distance; default 25"),
  make_option(c("-T","--threads"), default = 4, help = "number of threads")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

qFile <- args$qFile
tFile <- args$tFile
sFile <- args$sFile
oFile <- args$oFile

MIN_FRAC <- MAX_FRAC <- NULL

MIN_LEN <- as.numeric(args$min_len)
if(MIN_LEN < 1){
  MIN_FRAC <- MIN_LEN
}
MAX_LEN <- as.numeric(args$max_len)
if(MAX_LEN < 10){
  MAX_FRAC <- MAX_LEN
}


print(paste("Using min length", MIN_LEN, "; max length", MAX_LEN, "; min fraction", MIN_FRAC, "; max fraction", MAX_FRAC))

PROTECT_DIST <- args$protect_dist
# testing -----------------------------------------------------------------

# setwd("/mnt/mdwilson/huayun/test/test_halper/")
# 
# qFile <- "maec_rela_peaks.bed" #original peak file
# tFile <- "mouse_to_human_halLiftover_out.bed" # halLiftover output bed file
# sFile <- "mouse_to_human_halLiftover_summit.out.bed"   # summits halLiftover output file
# oFile <- "mouse_to_human_halper_Min100Max1500Protect25.bed" # output file
# 
# MIN_LEN <- 100
# MAX_LEN <- 1000
# PROTECT_DIST <- 25

# read in files
summits <- read.table(sFile, sep = "\t", col.names = c("seqnames","start","end","summit_peakname"))
targets <- read.table(tFile, sep = "\t", col.names = c("seqnames","start","end","mapped_peakname"))
peaks <- read.table(qFile, sep = "\t", col.names = c("seqnames","start","end","peakname"))

peak_width <- peaks %>% 
  mutate(width = end - start) %>% 
  dplyr::select(peakname, width) %>% 
  deframe()

summits_gr <- makeGRangesFromDataFrame(summits, keep.extra.columns = T)
targets_gr <- makeGRangesFromDataFrame(targets, keep.extra.columns = T)
targets_gr_list <- split(targets_gr, factor(targets_gr$mapped_peakname))



# overlap summits and targets -- not necessary
# ol <- findOverlaps(summits_gr, targets_gr_list)
# ol_peaknames <- cbind(summits[from(ol),],
#                       mapped_peakname = names(targets_gr_list)[to(ol)]) %>% 
#   filter(summit_peakname == mapped_peakname)

# every mapped summit will correspond to a mapped region by default
# go through summit regions one by one

chain_ortho_region <- function(summit_used){
 # print(summit_used)
  width <- peak_width[[summit_used$summit_peakname]]
  
  if(!is.null(MAX_FRAC)){
    max_len <- width * MAX_FRAC
  } else{
    max_len <- MAX_LEN
  }
  
  if(!is.null(MIN_FRAC)){
    min_len <- width * MIN_FRAC
  } else{
    min_len <- MIN_LEN
  }
  
  summit_used <- makeGRangesFromDataFrame(summit_used, keep.extra.columns = T)
  summit_extended <- resize(summit_used, width = 2*max_len, fix = "center")
  
  fragments <- GenomicRanges::reduce(targets_gr_list[[summit_used$summit_peakname]]) # merge adjacent fragments
  fragments <- sort(subsetByOverlaps(fragments, summit_extended)) # filter for only fragments that are on the same chromosome as the mapped summit, always sort by chromosome coordinates
  
  ol <- findOverlaps(fragments, summit_used)  # find the fragment the summit maps to
  start_index <- from(ol) # starting from the fragment that summit maps to
  
  get_length <- function(left, right){
    # given fragments, left and right index, calculate length of chained regions
    new_start <- start(fragments[left])
    new_end <- end(fragments[right])
    width <- new_end - new_start
    summit_start <- start(summit_used) - new_start
    summit_right <- new_end - start(summit_used)
    return(c("chr" = as.character(seqnames(summit_used)), "start" = new_start, "end" = new_end, "summit" = start(summit_used), "peakname" = summit_used$summit_peakname, "width" = width, "summit_start_dist" = summit_start, "summit_end_dist" = summit_right))
  }
  
  if(length(fragments) == 1){
    final <- as.data.frame(t(get_length(1,1))) %>% 
      mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len & as.numeric(summit_start_dist) >= PROTECT_DIST & as.numeric(summit_end_dist) >= PROTECT_DIST)
  } else{
    all_combs <- as.data.frame(t(combn(1:length(fragments), 2))) %>% 
      setNames(c("left","right")) %>% 
      filter(left <= start_index & right >= start_index) %>% 
      add_row(left = start_index, right = start_index)
    
    out_matrix <- as.data.frame(t(apply(all_combs, 1, function(x) get_length(x[1],x[2])))) %>% 
      mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len & as.numeric(summit_start_dist) >= PROTECT_DIST & as.numeric(summit_end_dist) >= PROTECT_DIST)
    
    # pick the best match, 
    final <- out_matrix %>% 
      filter(filter) %>% 
      filter(width == max(width))
    if(nrow(final) == 0){
      final <- out_matrix %>% 
        filter(width == max(width))
    }
  }
  final$ori_width <- width
  return(final)
}


ortho_regions <- do.call(rbind, mclapply(1:nrow(summits), function(x) chain_ortho_region(summits[x,]), mc.cores = args$threads))

#saveRDS(ortho_regions, paste0(oFile, ".rds"))

#ortho_regions <- do.call(rbind, mclapply(1:100, function(x) chain_ortho_region(summits[x,]), mc.cores = 2))

#ortho_regions <- lapply(1:nrow(summits), function(x) chain_ortho_region(summits[x,]))
ortho_regions <- dplyr::select(ortho_regions, chr, start, end, summit, peakname, ori_width, width, summit_start_dist, summit_end_dist, filter)

write.table(subset(ortho_regions, filter == TRUE) %>% dplyr::select(-filter), file = oFile, sep = "\t", row.names = F, col.names = F, quote = F)
write.table(subset(ortho_regions, filter == FALSE) %>% dplyr::select(-filter), file = paste0(oFile,".failed"), sep = "\t", row.names = F, col.names = F, quote = F)
