#!/usr/bin/Rscript

#This script combines output from summary statistics generated by Samtools and Bedtools
#It requires input in tab-separated format 


#Import libraries
library(methods)

#Set up argument parser and QC-------------------------------------------------
args <- commandArgs()

#print(args)
cov <- args[6]
read <-args[7]
output <-args[8]

if(is.na(cov)) 
{
  stop("Coverage file not defined.")
}

if(is.na(read)) 
{
  stop("Mapped read file not defined.")
}

if(is.na(output)) 
{
  stop("Output file not defined")
}

#Define sample name------------------------------------------------------------
# cov <- "test-A.genomecoverage.txt"
# read <- "test-A.sorted.idxstats.tsv"
sample_name = sub(".genomecoverage.txt", "", cov)
sample_name = sub(".*\\/", "", sample_name)

##Build a dataframe of coverage and mapped reads-------------------------------
cov_df <- read.delim(cov, sep = "\t", header=F)
cov_df$V3 <- sample_name
colnames(cov_df) <- c("AlignTarget", "Coverage", "Sample")

read_df <-read.delim(read, sep = "\t", header=F)
read_df$V5 <-sample_name
colnames(read_df) <- c("AlignTarget", "TargetLength", 
                       "MappedReads", "UnmappedReads", 
                      "Sample")
read_df <- subset(read_df, AlignTarget != "*", 
                  select = c(Sample, AlignTarget, TargetLength, MappedReads))

#Merge data. Keep all observations
all_align_data <- merge(cov_df, read_df, 
                        by = c("Sample", "AlignTarget"), 
                        all.x = TRUE, all.y = TRUE)

##Write output
write.table(all_align_data, sep = "\t",
            file = output,
            quote = FALSE,
            col.names = FALSE, row.names = FALSE)

