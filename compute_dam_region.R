#!/usr/bin/env Rscript

# Load required libraries
library(argparse)
suppressMessages(library(data.table))
suppressMessages(library(rtracklayer))

# Create an argument parser
parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

# Add arguments to the parser
parser$add_argument('--counts', help='count file with counts per region')
parser$add_argument('--bed', help='bed file with regions of interest')
parser$add_argument('--window', help='window around integration')
parser$add_argument('--out', help='output')
parser$add_argument('--threads', help='number of threads', type="integer")
parser$add_argument('--bw', help='bigwig file with binned genome track')

# Parse the command line arguments
argv <- parser$parse_args()

# Set the number of threads
setDTthreads(threads=argv$threads)

# Save the parsed arguments to a file
save(argv, file='test.Rdata')

# Split the count file path into a vector
count_vec = unlist(strsplit(argv$counts, ','))

# Read the count files into a list of data tables
count_list = lapply(count_vec, function(x){
                        fread(cmd=paste('zcat', x), stringsAsFactors=F,
                              col.names=c('seqnames', 'start', 'end', 'count'))})

# Calculate the total counts per chromosome
chr_count = lapply(count_list, function(dt){dt[,.(total=sum(count)),by=.(seqnames)]})

# Import the bed file (integration site) as a genomic ranges object
map_gr = import.bed(argv$bed)

# Create a genomic ranges object for the GATC regions
gatc_gr = GRanges(count_list[[1]]$seqnames, IRanges(count_list[[1]]$start,
                                                    count_list[[1]]$end))

# Find overlaps between the GATC regions and the integration sites
o = findOverlaps(map_gr, gatc_gr, maxgap=as.numeric(argv$window)/2)

# Extract the GATC regions that overlap with the integration sites
gatc_list = lapply(count_list, function(x){
    dt = data.table(name=map_gr$name[queryHits(o)],
                    seqnames=x$seqnames[subjectHits(o)],
                    x[subjectHits(o), 'count'])
    setkey(dt, name)
    sum_dt = dt[,.(seqnames=seqnames[1], sum=sum(count)), by=.(name)]
})

# Calculate the total counts per chromosome for the GATC regions
total_list = lapply(1:length(gatc_list), function(i){
    dt = chr_count[[i]][gatc_list[[i]], on=.(seqnames)]
    dt$total
})

# Create a data table with the GATC region names and their counts
gatc_dt = data.table(name=gatc_list[[1]]$name,
                     do.call(cbind, lapply(gatc_list, function(x){x$sum})),
                     do.call(cbind, total_list))

# Define a function to calculate the geometric mean
gm_mean = function(x){
      prod(x)^(1/length(x))
}

# Calculate the normalized counts for each GATC region
# This script computes the normalized values, mean values, and z-scores for a given dataset.

# Calculate the number of regions
R = length(count_vec)/2

# Normalize the data
norm_list = lapply(1:R, function(i){
    i_vec = 0:3 * R + 1 + i
    ct = gatc_dt[,.SD,.SDcols=c(1,i_vec)]
    colnames(ct)[2:5] = c('exp', 'ctrl', 't_exp', 't_ctrl')
    norm_c=ct[,.(exp=exp/t_exp * min(t_exp, t_ctrl) + 1,
                 ctrl=ctrl/t_ctrl * min(t_exp, t_ctrl) + 1), by=name]
    norm_c[,list(norm=exp/ctrl, ctrl=ctrl)]
})

# Create a table of normalized values
norm_table = data.table(name=gatc_dt$name,
                        do.call(cbind, lapply(norm_list, function(x){x[,'norm']})))

# Calculate the mean and geometric mean for each region
mean_table = norm_table[,.(mean=mean(log2(unlist(.SD))),
                           gm_mean = gm_mean(log2(.SD))), by=name]

# Create a table of control values
ctrl_table = data.table(name=gatc_dt$name,
                        do.call(cbind, lapply(norm_list, function(x){x[,'ctrl']})))

# Calculate the mean and geometric mean for control values
mean_ctrl_table = ctrl_table[,.(mean_ctrl = mean(log2(unlist(.SD))),
                                gm_mean_ctrl = gm_mean(log2(.SD))), by=name]

# Merge the mean and control tables
mean_table = merge(mean_table, mean_ctrl_table, by='name')

# Import the track data
track_gr = import.bw(argv$bw)

# Calculate the mean and standard deviation of the track scores
mean_track = mean(track_gr$score)
sd_track = sd(track_gr$score)

# Calculate the z-score for each region
mean_table$z_score = (mean_table$mean - mean_track) / sd_track

# Write the mean table to a file
write.table(mean_table, file=argv$out, quote=F, sep='\t', row.names=F)
