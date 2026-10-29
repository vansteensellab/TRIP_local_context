# TRIP local context pipeline

# Table of Contents

- [**1. Introduction**](#1-introduction)
- [**2. Requirements**](#2-requirements)
- [**3. Input formats**](#3-input-formats)
  - [3.1 Integration sites (summary file)](#31-integration-sites-summary-file)
  - [3.2 Coverage tracks](#32-coverage-tracks)
  - [3.3 pA-DamID count files](#33-pa-damid-count-files)
  - [3.4 Feature files](#34-feature-files)
- [**4. Main steps of the pipeline**](#4-main-steps-of-the-pipeline)
  - [4.1 Expanding integration sites into windows](#41-expanding-integration-sites-into-windows)
  - [4.2 Coverage](#42-coverage)
  - [4.3 pA-DamID](#43-pa-damid)
  - [4.4 Overlapping and nearest features](#44-overlapping-and-nearest-features)
- [**5. Setting up a config file for the pipeline**](#5-setting-up-a-config-file-for-the-pipeline)
  - [5.1 Integration sites and windows](#51-integration-sites-and-windows)
  - [5.2 Data sources](#52-data-sources)
  - [5.3 Other required parameters](#53-other-required-parameters)
- [**6. Output files**](#6-output-files)
- [**7. Running the pipeline**](#7-running-the-pipeline)

# 1. Introduction

TRIP (Thousands of Reporters Integrated in Parallel) measures the activity of a reporter at thousands of different genomic locations. To interpret those measurements, you usually want to know what the chromatin *around* each integration site looks like.

This repo contains a Snakemake pipeline that annotates the **local genomic context** of every integration site of a TRIP experiment. For each integration site, plus a window of configurable size around it, the pipeline reports three types of information:

1. **Coverage**: the mean signal of any bigWig track (e.g. ChIP-seq, DNase-seq) in the window.
2. **pA-DamID signal**: the log2 ratio between a pA-DamID experiment and its control in the window, together with a z-score relative to the genome-wide distribution.
3. **Genomic features**: which features (genes, enhancers, etc.) overlap the window, and which feature is the nearest to the integration site.

# 2. Requirements

This repository does not ship a conda environment, so the following tools need to be available in your `PATH` (e.g., in a conda environment) before running the pipeline.

- **Snakemake**
- **deepTools** (`multiBigwigSummary`)
- **bedtools**
- **UCSC tools** (`bedGraphToBigWig`)
- Standard Unix tools: `awk`, `sort`, `zcat`
- **R** with the packages: `argparse`, `data.table`, `dplyr`, `GenomicRanges`, `rtracklayer`

# 3. Input formats

The pipeline requires an **integration sites file** and a **chromosome sizes file**. Everything else is optional: provide any combination of coverage tracks, pA-DamID experiments and feature files, depending on which information you want.

## 3.1 Integration sites (summary file)

A tab-separated file **with header**, one row per integration. It must contain the barcode, the chromosome, and the position of the integration:

| barcode | chr  | start    |
|---------|------|---------:|
| ACGT... | chr1 | 3014522  |
| TTGA... | chr7 | 98122340 |

**Column descriptions:**

- **barcode**: barcode of the integration. This column **must be named `barcode`**, as the scripts use this name to retrieve it.
- **chr**: chromosome of the integration (`chr`, `chrom` or `chromosome` are recognized).
- **start**: position of the integration.

Integrations on `chrM` are dropped.

## 3.2 Coverage tracks

Any number of bigWig files. Each one is summarized (mean signal) in every window around every integration site.

## 3.3 pA-DamID count files

For each pA-DamID experiment you need to provide the experiment (Dam-fused protein, e.g. Dam-H3K9me3) and control (e.g. Dam-only) files.
Each file is a gzipped, tab-separated file with the count per GATC fragment, with no header and four columns:

| chr  | start | end  | count |
|------|------:|-----:|------:|
| chr1 | 3000  | 3200 | 12    |

> The order of the files in the config matters: replicate 1 of the experiment is matched with replicate 1 of the control, and the same for replicate 2.

## 3.4 Feature files

Any number of `bed`, `gff` or `vcf` files with the genomic features you want to relate the integrations to (genes, enhancers, etc.).

# 4. Main steps of the pipeline

The pipeline is run separately for each site (a named set of integrations) and each window size. 
The three categories of information (coverage, pA-DamID, features) are independent of each other.

## 4.1 Expanding integration sites into windows

The integration positions from the summary file are turned into BED files, in which each integration is expanded to a window of the desired size, centered on the integration position. Windows are trimmed so that they do not exceed the chromosome sizes, and the output is sorted with `bedtools sort`.

- One BED file is created for each window size listed under `bins` (used for coverage and pA-DamID).
- One BED file is created with the window defined by `bin_for_overlap_feature` (used for overlapping features).
- One BED file is created with the window defined by `bin_for_nearest_feature` (used for nearest features). By default this is 1 bp, so that the exact position of the integration is used.

## 4.2 Coverage

For every window size, `multiBigwigSummary BED-file` computes the mean signal of every bigWig listed under `data_sources: coverage` in every window. The result is a single table with one row per window and one column per track.

## 4.3 pA-DamID

1. **Build tracks.** The GATC count files of each replicate (experiment and control) are converted into bedGraph and then into bigWig files.
2. **Bin tracks.** The four bigWig files (2 experiment replicates and 2 control replicates) of each experiment are binned into the genome-wide bin sizes defined under `bins`, using `multiBigwigSummary bins`.
3. **Normalize tracks.** For every bin and every replicate, the experiment and the control are scaled to the same total (the smaller of the two totals), a pseudocount of 1 is added, and the ratio experiment/control is computed. The final track score is the average of the log2 ratios across the two replicates. This is exported as a bigWig file.
4. **Get means.** For every integration site, the script sums the GATC counts in the window around the integration, normalizes experiment and control in the same way as above (using the total counts of the chromosome), and reports the log2 ratio. The z-score is calculated by comparing the value of each integration to the mean and standard deviation of the genome-wide normalized track from step 3.

## 4.4 Overlapping and nearest features

For every feature file listed under `data_sources: feature`:

- `bedtools intersect -wo` reports which features overlap the window defined by `bin_for_overlap_feature`, and the number of overlapping bases.
- `bedtools closest -d` reports the nearest feature to the window defined by `bin_for_nearest_feature` and its distance.

# 5. Setting up a config file for the pipeline

As in any Snakemake workflow, you first need to set up a configuration file. A full example is available in `config.yaml`. Below we guide you through the fields of the config file.

## 5.1 Integration sites and windows

```yaml
# Directory where all the results will be saved
output_directory: 'my_output_directory'

# Integration sites. Use one name per set of integrations (we call it 'site')
summary_file:
  my_TRIP: '/path/to/integration_coordinates.txt'

# Window sizes (in bp) to be used for each site.
# The names must match the ones used in 'summary_file'
bins:
  my_TRIP: [2000, 5000, 10000, 20000, 50000, 100000]

# Window (in bp) used to find overlapping features
bin_for_overlap_feature: 1000000

# Window (in bp) used to find the nearest feature.
# 1 means the exact integration position
bin_for_nearest_feature: 1
```

> You can add more than one site under `summary_file` and `bins`. The pipeline will produce one set of results per site.

## 5.2 Data sources

All three sub-sections are used to define the data that will be linked to the integrations. They can be given in any number. Names given here (e.g. `H3K9me3`) are used as column labels or file names in the output.

```yaml
data_sources:
  # pA-DamID experiments. Each one needs exactly two replicates of experiment and control
  pADamID:
    H3K9me3:                    # Name of the experiment
      experiment:
        - /path/to/rep1_Dam_k9me3-gatc.counts.txt.gz
        - /path/to/rep2_Dam_k9me3-gatc.counts.txt.gz
      control:
        - /path/to/rep1_Dam_only-gatc.counts.txt.gz
        - /path/to/rep2_Dam_only-gatc.counts.txt.gz

  # bigWig files. Name: path
  coverage:
    H3K9me3: /path/to/H3K9me3.bigWig
    H3K27me3: /path/to/H3K27me3.bigWig
    DNase: /path/to/DNase.bigWig

  # Feature files (bed, gff or vcf). Name: path
  feature:
    gene: /path/to/genes.gff3
    enhancers: /path/to/enhancers.bed
```

## 5.3 Other required parameters

```yaml
# Chromosome sizes file (two columns: chromosome and length)
chrom_sizes:
  mm10: /path/to/mm10.chrom.sizes
```

# 6. Output files

```
output_directory/
├── coverage_results/        # <- Mean signal of every coverage track in each window
│   └── {site}-{binsize}_summary.txt
├── pADamID_results/         # <- Final pA-DamID results per integration
│   └── {site}-{binsize}_{experiment}.txt
├── feature_results/         # <- Overlapping and nearest features
│   ├── overlap_{site}_{feature}.txt
│   └── nearest_{site}_{feature}.txt
├── inputs/                  # <- Windows around integrations (BED)
├── feature_input/           # <- Windows used for overlap and nearest features (BED)
├── pADamID_tracks/          # <- Intermediate: per-replicate bigWig tracks
├── pADamID_bins/            # <- Intermediate: binned counts
└── pADamID_normalized_tracks/   # <- Intermediate: normalized log2 ratio tracks (bigWig)
```

**`coverage_results/*_summary.txt`**: one row per window (`chr`, `start`, `end`) and one column per coverage track, with the mean signal in that window.

**`pADamID_results/*.txt`**: one row per integration (`name` = barcode), with the columns:

| Column | Description |
|---|---|
| `mean` | Mean across replicates of the log2(experiment / control) |
| `gm_mean` | Geometric mean of the same per-replicate values |
| `mean_ctrl` | Mean across replicates of the log2 normalized control |
| `gm_mean_ctrl` | Geometric mean of the same control values |
| `z_score` | `mean` standardized against the genome-wide distribution of the normalized track |

**`feature_results/overlap_*.txt`**: output of `bedtools intersect -wo`. The first four columns are the window (chr, start, end, barcode), followed by the columns of the overlapping feature, and the number of overlapping bases in the last column.

**`feature_results/nearest_*.txt`**: output of `bedtools closest -d`. The first four columns are the integration (chr, start, end, barcode), followed by the columns of the nearest feature, and the distance in bp in the last column.

# 7. Running the pipeline

1. **Set up your configuration file**: Copy and modify `config.yaml` with your parameters.
2. **Activate your environment**: Make sure all the tools listed in [Requirements](#2-requirements) are available.
3. **Dry run** (recommended): Test your configuration without running the pipeline:

   ```bash
   snakemake --configfile your_config.yaml --dry-run
   ```
4. **Run the pipeline**: Execute from your analysis folder with your desired number of cores:

   ```bash
   snakemake --cores 80 --printshellcmds --configfile your_config.yaml
   ```

If you cloned this repository elsewhere, point Snakemake to the Snakefile with `-s /path/to/TRIP_local_context/Snakefile`. 
The R scripts are found automatically, as they are expected to be in the same folder as the Snakefile.
