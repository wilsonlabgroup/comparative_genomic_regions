#!/usr/bin/env Rscript
# implements HALPER method for extending hal mapped summits

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
summits <- read.table(sFile, sep = "\t")[,1:4] %>% setNames(c("seqnames","start","end","summit_peakname")) # lifted over summits
targets <- read.table(tFile, sep = "\t", col.names = c("seqnames","start","end","mapped_peakname")) # lifted over peaks
peaks <- read.table(qFile, sep = "\t", col.names = c("seqnames","start","end","peakname")) # original peaks

# get the width of the original peaks
peak_width <- peaks %>% 
  mutate(width = end - start) %>% 
  dplyr::select(peakname, width) %>% 
  deframe()

# convert regions to granges 
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

chain_ortho_region <- function(peak_used, summit = TRUE){
  # if its summit, input 'peak_used' is a line from mapped summit file; else input is a line from mapped peaks 
  if(summit){
    peakname <- peak_used$summit_peakname
  } else{
    #peakname <- unlist(strsplit(peak_used,"-"))[1]
    peakname <- peak_used$mapped_peakname[1]
  }
  width <- peak_width[[peakname]]
  
  print(peakname)
  
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
  
  # get a grange of fragments to be chained 
  if(summit){
    summit_used <- makeGRangesFromDataFrame(peak_used, keep.extra.columns = T)
    summit_extended <- resize(summit_used, width = 2*max_len, fix = "center")
    
    fragments <- GenomicRanges::reduce(targets_gr_list[[peak_used$summit_peakname]]) # merge adjacent fragments
    fragments <- sort(subsetByOverlaps(fragments, summit_extended)) # filter for only fragments that are on the same chromosome as the mapped summit, always sort by chromosome coordinates
    
    ol <- findOverlaps(fragments, summit_used)  # find the fragment the summit maps to
    start_index <- from(ol) # starting from the fragment that summit maps to
  } else{
      # input is multiple fragments mapped from a peak
    
      fragments <- GenomicRanges::reduce(peak_used)
    
      start_index <- which.max(width(fragments)) # get longest fragment and use it as an anchor
      longest_extended <- resize(fragments[start_index], width = 2*max_len, fix = "center")
      fragments <- sort(subsetByOverlaps(fragments, longest_extended))
      start_index <- which.max(width(fragments)) # get index for longest fragment again
  }
  
  
  get_length <- function(left, right){
    # given fragments, left and right index, calculate length of chained regions
    new_start <- start(fragments[left])
    new_end <- end(fragments[right])
    width <- new_end - new_start
    if(summit){
      summit_start <- start(summit_used) - new_start
      summit_right <- new_end - start(summit_used)
      return(c("chr" = as.character(seqnames(summit_used)), "start" = new_start, "end" = new_end,  "peakname" = peakname,"width" = width, "summit" = start(summit_used),  "summit_start_dist" = summit_start, "summit_end_dist" = summit_right))
    } else{
      return(c("chr" = as.character(seqnames(fragments[1])), "start" = new_start, "end" = new_end, "peakname" = peakname, "width" = width))
    }
  }
  

  if(length(fragments) == 1){
    if(summit){
      final <- as.data.frame(t(get_length(1,1))) %>% 
        mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len & as.numeric(summit_start_dist) >= PROTECT_DIST & as.numeric(summit_end_dist) >= PROTECT_DIST)
    } else{
      final <- as.data.frame(t(get_length(1,1))) %>% 
        mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len)
    }
  } else{
    all_combs <- as.data.frame(t(combn(1:length(fragments), 2))) %>% 
      setNames(c("left","right")) %>% 
      filter(left <= start_index & right >= start_index) %>% 
      add_row(left = start_index, right = start_index)
    
    if(summit){
      out_matrix <- as.data.frame(t(apply(all_combs, 1, function(x) get_length(x[1],x[2])))) %>% 
        mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len & as.numeric(summit_start_dist) >= PROTECT_DIST & as.numeric(summit_end_dist) >= PROTECT_DIST)
    }else{
      out_matrix <- as.data.frame(t(apply(all_combs, 1, function(x) get_length(x[1],x[2])))) %>% 
        mutate(filter = as.numeric(width) >= min_len & as.numeric(width) <= max_len)  
    }
    
    # pick the best match, 
    final <- out_matrix %>% 
      filter(filter) %>% 
      filter(width == max(width))
    if(nrow(final) == 0){
      final <- out_matrix %>% 
        filter(width == max(width))
    }
  }

  if(!"summit" %in% colnames(final)){
    final[, c("summit","summit_start_dist","summit_end_dist")] <- "."
  }
  final$ori_width <- width
  final <- final %>% 
    dplyr::select(chr, start, end, summit, peakname, ori_width, width, summit_start_dist, summit_end_dist, filter)
  if(final$filter){
    # if the chained region passed cutoff, return all actual fragments for overlapping with peaks
    final_range <- GRanges(final$chr, IRanges(as.numeric(final$start), as.numeric(final$end)))
    out_fragments <- as.data.frame(subsetByOverlaps(fragments, final_range))[,1:3]
    out_fragments$peakname <- peakname
  } else{
    out_fragments <- NULL
  }
  
  return(list(final, out_fragments))
}

temp1 <- mclapply(1:nrow(summits), function(x) chain_ortho_region(summits[x,]), mc.cores = args$threads)
ortho_regions <- do.call(rbind, lapply(temp1, "[[",1))
ortho_fragments <- do.call(rbind, lapply(temp1, "[[",2))

targets_no_summits_gr <- makeGRangesFromDataFrame(subset(targets, !mapped_peakname %in% summits$summit_peakname) %>% unite(new_peakname, mapped_peakname, seqnames, sep = "-", remove = F), keep.extra.columns = T)

targets_no_summits_gr_list <- split(targets_no_summits_gr, factor(targets_no_summits_gr$new_peakname))

temp2 <- mclapply(names(targets_no_summits_gr_list), function(x) chain_ortho_region(targets_no_summits_gr_list[[x]], summit = F), mc.cores = args$threads)
ortho_regions_summit_notMapped <- do.call(rbind, lapply(temp2, "[[",1))
ortho_fragments_summit_notMapped <- do.call(rbind, lapply(temp2, "[[",2))
#save.image("test.RData")

#ortho_regions_summit_notMapped <- lapply(names(targets_no_summits_gr_list), function(x) chain_ortho_region(targets_no_summits_gr_list[[x]], summit = F))


#saveRDS(ortho_regions, paste0(oFile, ".rds"))

#ortho_regions <- do.call(rbind, mclapply(1:100, function(x) chain_ortho_region(summits[x,]), mc.cores = 2))

#ortho_regions <- lapply(1:nrow(summits), function(x) chain_ortho_region(summits[x,]))
#ortho_regions <- dplyr::select(ortho_regions, chr, start, end, summit, peakname, ori_width, width, summit_start_dist, summit_end_dist, filter)

out <- rbind(subset(ortho_regions, filter == TRUE) %>% dplyr::select(-filter),
             subset(ortho_regions_summit_notMapped, filter == TRUE) %>% dplyr::select(-filter))
failed <- rbind(subset(ortho_regions, filter == FALSE) %>% dplyr::select(-filter),
                subset(ortho_regions_summit_notMapped, filter == FALSE) %>% dplyr::select(-filter))
out_fragments <- rbind(ortho_fragments, ortho_fragments_summit_notMapped)

suffix <- unlist(strsplit(oFile, "\\."))[-1]
write.table(out_fragments, gsub(paste0(".", suffix),"_fragments.bed",oFile), sep = "\t", col.names = F, row.names = F, quote = F) # output overlapping fragments

write.table(out, file = oFile, sep = "\t", row.names = F, col.names = F, quote = F) # output chained regions

write.table(failed, file = paste0(oFile,".failed"), sep = "\t", row.names = F, col.names = F, quote = F) # output failed regions 
