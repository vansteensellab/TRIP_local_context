#!/usr/bin/env Rscript

# Expand integration position to desired window size
# ==================================================
# Vinícius H. Franceschini Santos, 2024-06-12
#
# This script reads the TRIP summary file (a BED-like file) and reads the desire
# window size to expand the integration position. It rund Granges functions to
# perform this expansion. The output is a BED file with the expanded integration
# position. It also gets the chromosome sizes information, so that the expanded
# integration position does not exceed the chromosome size.

suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(dplyr))

parser <- ArgumentParser()
parser$add_argument(
    "--input_bed",
    help = "Input file (TRIP summary file)"
)
parser$add_argument(
    "--bin_size",
    help = "Window size to expand the integration position"
)
parser$add_argument(
    "--output_bed",
    help = "Output file (expanded TRIP summary file)"
)
parser$add_argument(
    "--genome",
    help = "Name of the genome (e.g. mm10)"
)
args <- parser$parse_args()

# Read input bed; columns are name, chrom, start
bed <- fread(args$input_bed, header=T)
bed$end = bed$start
genome = Seqinfo(genome = args$genome)


# Create the GRanges object
expanded_gr =
    makeGRangesFromDataFrame(
        bed,
        seqinfo = genome,
        keep.extra.columns = T
    ) %>%
    # Expand to the desired window size
    resize(
        width = as.numeric(args$bin_size),
        fix = "center"
    )

# Removing sequences from chrM -- avoiding future error
expanded_gr = expanded_gr[seqnames(expanded_gr) != "chrM"]

# Trim to the extent of the chromosomes
expanded_gr = trim(expanded_gr)

# Convert to data frame
expanded_gr %>%
    as.data.frame() %>%
    # Select only the columns we want (BED format)
    select(seqnames, start, end, barcode) %>%
    # Save the results
    write.table(
        file = args$output_bed,
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
    )
