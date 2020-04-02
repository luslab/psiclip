#!/usr/bin/env Rscript

# Charlotte Capitanchik 2018 #
# Get cDNA counts from bed file #

library(data.table)
library(dplyr)
options(scipen = 999)
options("scipen"=100, "digits"=4)
args = commandArgs(trailingOnly=TRUE) # get command line arguments in a list [args]

# Run command: Rscript --vanilla scripts/getcDNAcounts.R {input.bed} {output.cdnacounts}

# INPUT #
input_file <- args[1]
output_file <- args[2]

input_bed <- fread(input_file) # Read input file
input_bed$V4 <- gsub(".*rbc:","",input_bed$V4) # Remove the unique read identifier and leave the random barcode
  
input_bed_forward <- input_bed %>% filter(V6=="+")  # Filter based on strand
if (nrow(input_bed_forward)==0){
  input_bed_forward <- data.frame(as.character(c("EMPTY")),as.numeric(c(1)),as.numeric(c(1)),as.character(c(".")),as.numeric(c(1)), stringsAsFactors = FALSE)
  names(input_bed_forward) <- c("chr","start","stop","strand","count")
} else {
input_bed_forward <- input_bed_forward %>% group_by(V1,V2,V4,V6) %>% filter(row_number()==1) # Get the unique cDNAs, FORWARD STRAND, so using V2
input_bed_forward$V3 <- input_bed_forward$V2 # Make the coordinate the xlsite (xlsite; CDNA START - 1)
input_bed_forward$V2 <- input_bed_forward$V2 -1
input_bed_forward <- ungroup(input_bed_forward)
input_bed_forward$number <- 1
input_bed_forward <- input_bed_forward %>% group_by(V1,V2,V3,V6) %>% summarise(total=as.integer(sum(number))) # Sum the xlsites to get summed cDNA counts
names(input_bed_forward) <- c("chr","start","stop","strand","count")
input_bed_forward <- as.data.frame(input_bed_forward)
}

input_bed_reverse <- input_bed %>% filter(V6=="-")  # Filter based on strand
if (nrow(input_bed_reverse)==0){
  input_bed_reverse <- data.frame(as.character(c("EMPTY")),as.numeric(c(1)),as.numeric(c(1)),as.character(c(".")),as.numeric(c(1)), stringsAsFactors = FALSE)
  names(input_bed_reverse) <- c("chr","start","stop","strand","count")
} else {
input_bed_reverse <- input_bed_reverse %>% group_by(V1,V3,V4,V6) %>% filter(row_number()==1) # Get the unique cDNAs, REVERSE STRAND, so using V3
input_bed_reverse$V2 <- input_bed_reverse$V3
input_bed_reverse$V3 <- input_bed_reverse$V3 + 1  # Make the coordinate the xlsite (xlsite; CDNA START - 1)
input_bed_reverse <- ungroup(input_bed_reverse)
input_bed_reverse$number <- 1
input_bed_reverse <- input_bed_reverse %>% group_by(V1,V2,V3,V6) %>% summarise(total=as.integer(sum(number))) # Sum the xlsites to get summed cDNA counts
names(input_bed_reverse) <- c("chr","start","stop","strand","count")
input_bed_reverse <- as.data.frame(input_bed_reverse)
}

# OUTPUT #
# Finally merge the dataframes and print to file
final_bed <- rbind(input_bed_forward, input_bed_reverse) # merge
final_bed <- final_bed %>% filter(chr != "EMPTY")
final_bed <- data.frame(final_bed$chr, final_bed$start, final_bed$stop, ".", final_bed$count, final_bed$strand) # format as a bed file
write.table(final_bed, file=output_file, col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)
