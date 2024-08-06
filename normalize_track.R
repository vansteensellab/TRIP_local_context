#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(rtracklayer))

# Define command line arguments
parser <- ArgumentParser(description='Process counts into normalized means for a track')
parser$add_argument('--input', help='counts of individual experiments')
parser$add_argument('--chrom_sizes', help='chrom sizes')
parser$add_argument('--out', help='output')
parser$add_argument('--binsize', help='binsize used')
parser$add_argument('--threads', help='number of threads', type="integer")

# Parse command line arguments
argv <- parser$parse_args()

# Set number of threads
setDTthreads(threads=argv$threads)

# Save command line arguments to a file
save(argv, file='test.Rdata')

# Read input data
input_dt = fread(argv$input, sep='\t')

# Rename column names by removing quotes
colnames(input_dt) = gsub("[#]*'(.*)'", '\\1', colnames(input_dt))

# Remove rows with NA or sum of counts <= 0
na_sum = rowSums(input_dt[,-c('chr', 'start', 'end')])
input_dt = input_dt[!is.na(na_sum) & na_sum > 0, ]

# Convert input data to GRanges object
track_gr = makeGRangesFromDataFrame(input_dt)

# Multiply count data by binsize
count_dt = input_dt[,-c('chr', 'start', 'end')] * as.numeric(argv$binsize)

# Simple scaling normalization.
# This function normalizes 'exp' and 'ctrl' so that they have the 
# same total sum (equal to the minimum of the original sums), and then 
# it returns the ratio of 'exp' to 'ctrl'.
norm <- function(exp, ctrl){
    n_vec = c(sum(exp), sum(ctrl))
    m_min = min(n_vec)
    exp=exp/n_vec[1] * m_min + 1
    ctrl=ctrl/n_vec[2] * m_min + 1
    return(exp / ctrl)
}

# Apply normalization function to count data
norm_dt = count_dt[,list(r1=norm(exp1, ctrl1),
                         r2=norm(exp2, ctrl2))]

# Calculate scores based on normalized data
score(track_gr) = rowMeans(log2(norm_dt))

# Read chromosome sizes
seql_df = read.table(argv$chrom_sizes)

# Create sequence vector with chromosome sizes
seq_vec = seql_df[,2]
names(seq_vec) = seql_df[,1]

# Set sequence lengths in track_gr object
seqlengths(track_gr) = seq_vec[seqlevels(track_gr)]

# Adjust start positions by adding 1
start(track_gr) = start(track_gr) + 1

# Export track_gr as a bigWig file
export.bw(trim(track_gr), argv$out)
