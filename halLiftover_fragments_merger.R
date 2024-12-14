#!/usr/bin/env Rscript
# simple method to cancatenate halliftover outputs

# Merge hal liftover output fragments that are within certain bps (controled by --gap). If fragments are on multiple chromosomes treat each chromosome separately
# After merging, filter the merged fragmenets by length relative to original query peaks (controled by --max_frac and --min_frac)
# If a mapped summit file is given, further filter the merged regions for ones that overlapped corresponding mapped summit files  

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))

option_list <- list(
  make_option(c("-q","--qFile"), 
              help="Input peak file, first 4 columns must be in standard bed format, peak names must be unique"),
  make_option(c("-t","--tFile"), help = "Input mapped bed file"),
  make_option(c("-o","--oFile"), help = "Output file"),
  make_option(c("-s","--sFile"), help = "input mapped-summit bed file"),
  make_option("--max_frac",default = 1.25, help = "Maximum fraction of clustered regions/peaks width, Default 125 (%). If >10, treated as maximum length (bp)"),
  make_option("--min_frac", default = 0.75, help = "Minimum fraction of clustered regions/peaks width, Default 25 (%). If >10, treated as minumum length (bp)"),
  #make_option("--max_len", help = "Maximum length of clustered regions/peaks width"),
  #make_option("--min_len", default = 75, help = "Minimum length of clustered regions/peaks width"),
  make_option("--gap", default = 50, help = "Gap width. Fragments within distance < gap width will be merged. Default 50")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# get parameters
qFile <- args$qFile
tFile <- args$tFile
oFile <- args$oFile

MAX_GAP <- args$gap
MIN_FRAC <- args$min_frac
MAX_FRAC <- args$max_frac

print(args)

# if either min/max frac is > 10 , treat it as min/max length
if(MIN_FRAC > 10 | MAX_FRAC > 10){
  filter_by <- "length"
  MAX_LEN <- args$max_frac
  MIN_LEN <- args$min_frac
  print(paste0("Filter merged regions by fixed length, between ", MIN_LEN, " and ", MAX_LEN))
} else{
  filter_by <- "frac"
  print(paste0("Filter merged regions by relative fraction of original peaks, between ", MIN_FRAC, " and ", MAX_FRAC))
}  

# check if mapped summit files exist
if(!"sFile" %in% names(args)){
  print("No mapped summits given, ignore summit overlapping requirement")
  use_summit <- FALSE
} else{
  sFile <- args$sFile
  use_summit <- TRUE
  print(paste0("Summit file:", sFile, " given. Require summits overlapping"))
}

############################################################
# # for testing
# MAX_GAP = 50
# MIN_FRAC = 0.75
# MAX_FRAC = 1.25
#  
# setwd("~/mdwilson/huayun/test/test_hal_50bp//")
# qFile <- "../test_halper/maec_rela_peaks.bed" #original peak file
# tFile <- "../test_halper/mouse_to_human_halLiftover_out.bed" # halLiftover output bed file
# oFile <- "mouse_to_human_halper_gap50_perc25_125.bed"
############################################################


# read in files
# lifted over peaks. Add a column "peakname_new" which addes chr names to peak name for easy parsing 
targets <- read.table(tFile, sep = "\t", col.names = c("seqnames","start","end","mapped_peakname")) %>% mutate(mapped_peakname_new = paste(mapped_peakname, seqnames, sep = ".")) 

# original peaks
peaks <- read.table(qFile, sep = "\t", col.names = c("seqnames","start","end","peakname")) 


# get the width of the original peaks
peak_width <- peaks %>% 
  mutate(width = end - start) %>% 
  dplyr::select(peakname, width) %>% 
  deframe()

# convert hal liftover target fragments to a granges list. Each grange object in the list is mapped fragments for one peak on one chromosome
targets_gr_list <- makeGRangesListFromDataFrame(targets, split.field = "mapped_peakname_new", keep.extra.columns = T)

# merge fragments within N (default=50) bp
targets_gr_list_merged <- GenomicRanges::reduce(targets_gr_list, min.gapwidth = MAX_GAP)

# convert the merged gr list object to data frame
temp <- unlist(targets_gr_list_merged)
peaknames <- names(temp)
targets_merged_df <- as.data.frame(temp, row.names = NULL) 
targets_merged_df$mapped_peakname_new <- peaknames
targets_merged_df$mapped_peakname <- sapply(strsplit(peaknames, "\\."),"[[",1)

# associate with original peaks width
targets_merged_df$peak_width <- peak_width[targets_merged_df$mapped_peakname]

# calculate percentage of mapped region clusters / peak width
targets_merged_df$perc <- with(targets_merged_df, width/peak_width)

# if mapped summit file exist

if(use_summit){
  # read in summit 
  summits <- read.table(sFile, sep = "\t")[,1:4] %>% setNames(c("seqnames","start","end","summit_peakname")) %>% makeGRangesFromDataFrame(keep.extra.columns = T) # lifted over summits
  # overlap mapped summits and merge fragments
  ol <- findOverlaps(summits, targets_gr_list_merged)
  # get the summit-peak name overlap
  ol_names <- data.frame(summit_name = summits$summit_peakname[from(ol)],
                         mapped_peakname_new = names(targets_gr_list_merged[to(ol)]))
  # add mapped summits name to the merged fragments data frame and filter only matching ones
  targets_merged_df <- left_join(targets_merged_df, ol_names, relationship = "many-to-many") %>% 
    filter(mapped_peakname == summit_name)
}



# filter merged regions based on size (> merging gap threshold and within a certain proportion of the original peak) or by merged width

if(filter_by == "frac"){
  targets_merged_df_fil <- subset(targets_merged_df, perc >= MIN_FRAC & perc <= MAX_FRAC)
} else{
  targets_merged_df_fil <- subset(targets_merged_df, peak_width >= MIN_LEN & peak_width <= MAX_FRAC)
}


#### filter the original fragments for only the ones that are merged
# convert df back to grangeslist
targets_merged_fil_list <- makeGRangesListFromDataFrame(targets_merged_df_fil, split.field = "mapped_peakname_new")

# subset original fragments granges by peak names
ori_targets <- targets_gr_list[names(targets_merged_fil_list)]

# filter original fragments by overlapping with merged fragments
ori_targets_fil <- subsetByOverlaps(ori_targets, targets_merged_fil_list)


##### format for output
out <- targets_merged_df_fil[,c("seqnames","start","end","mapped_peakname")] %>% 
  arrange(mapped_peakname)
out_fragments <- as.data.frame(unlist(ori_targets_fil), row.names = NULL)[,c("seqnames","start","end","mapped_peakname")] %>% 
  arrange(mapped_peakname)

# wirte output 
suffix <- unlist(strsplit(oFile, "\\."))[-1]
if(length(suffix) == 0){
  out_frag_file <- paste0(oFile,"_fragments.bed")
} else{
  out_frag_file <- gsub(paste0(".", suffix),"_fragments.bed",oFile)
}
write.table(out_fragments, file = out_frag_file, sep = "\t", col.names = F, row.names = F, quote = F) # output overlapping fragments

write.table(out, file = oFile, sep = "\t", row.names = F, col.names = F, quote = F) # output chained regions





