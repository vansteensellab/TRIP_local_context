###############################################################################
## TRIP location info pipeline ################################################
###############################################################################
## Vinícius H. Franceschini-Santos
## 2024-06-11
## version 1.0

# Introduction snakemake:
# Snakemake is a bioinformatics pipeline tool written in python. It follows the
# 'make'-logic: you determine rules with which files can be generated, and
# which input files these require. Snakemake will try to execute the first rule
# given, and execute other rules based on the input the final rule requires.
# In this way, you have a pipeline.
# Some cool features of snakemake:
# - You give input and output file names, and snakemake will do the rest.
# - If a rule gives an error, the output file is automatically removed, to
#   prevent "partial" data files.
# - Given a number of cores to snakemake, it will automatically parallelize
#   things when possible.
# - You can create a fancy workflow overview picture.

# This is rework of a pipeline originally written by Christ. The original pipeline
# is in /DATA/scratch/usr/c.leemans/projects/DSB_repair_TRIP/pipeline/location_info/
# This pipeline is designed to provide some information about the location of
# the integration sites of a TRIP experiment. The three main categories of information
# are: (i) coverage of bigwig files, for any given file, (ii) difference in
# pA-DamID signal, for any given pA-DamID experiment, (iii) distance to the
# nearest featuer, for any given gene or other feature file.
# By 'location', we mean the integration site of the TRIP experiment plus a
# window around it. You can set different window sizes (bins) in the config file.

# Versions
# 1.0 - Initial version, vhfsantos, 2024-06-11

# Execute as (from analysis folder):
# snakemake --cores 80 --printshellcmds --configfile your_config_file.yaml


###############################################################################
### General set-up ############################################################
###############################################################################

import inspect
import os

filename = inspect.getframeinfo(inspect.currentframe()).filename
scripts_dir = os.path.dirname(os.path.abspath(filename))

##################################
### Helper functions #############
##################################
# Here I'll define some helper functions. These are mainly used to get the
# input files for the rules. Some of them also define the wildcards for the
# rules.


def get_coverage_file(config):
    """
    Get the bigwig files for the coverage tracks.
    ---
    Used in the coverage_summarize_results rule.
    """
    cov_dict = config["data_sources"]["coverage"]
    for key in cov_dict:
        yield (cov_dict[key])


def get_coverage_ID(config):
    """
    Get the IDs for the coverage tracks.
    ---
    Used in the coverage_summarize_results rule.
    """
    cov_dict = config["data_sources"]["coverage"]
    for key in cov_dict:
        yield (key)


def get_beds(config, wildcards):
    """
    Returns the bed file for the given site and bin size.
    Defines the wildcards:
    * site (the title of the summary file, defined in the config)
    * binsize (the windows defined in the config: 2000, 5000, 10000, ...).
    ---
    Used in the rules: coverage_summarize_results, paDamID_get_means.
    """
    return "%s/inputs/trimmed_%s-%s.bed" % (
        wildcards.outdir,
        wildcards.site,
        wildcards.binsize,
    )


def get_paDamID_files(config, wildcards):
    """
    This function returns the pA-DamID files for the given experiment and condition.
    Defines the wildcards:
    * dam_experiment (title of  the pA-DamID experiment, defined in the config: Lamin, H3K9me3,.. )
    * condition (condition of each file: experiment or control)
    * rep (replicate number: 1, 2, 3, ...)
    ---
    Used in the paDamID_build_bigwig_track rule.
    """
    # Subtract 1 because list indexing starts at 0, but replicate numbers start at 1.
    # this makes the 'rep' as wildcard 1-based
    index = int(wildcards.rep) - 1
    return config["data_sources"]["pADamID"][wildcards.dam_experiment][
        wildcards.condition
    ][index]


def get_paDamID_counts(config, wildcards):
    """
    Returns the count files for the given pA-DamID experiment and condition.
    Defines the wildcards:
    * dam_experiment (title of the pA-DamID experiment, defined in the config: Lamin, H3K9me3,.. )
    ---
    Used in the paDamID_build_bigwig_track rule.
    """
    dam_dict = config["data_sources"]["pADamID"][wildcards.dam_experiment]
    return [exp for exp in dam_dict["experiment"]] + [
        dam for dam in dam_dict["control"]
    ]


def get_paDamID_normalized_tracks(config, wildcards):
    """
    Returns the track file for the given pA-DamID experiment and condition.
    Defines the wildcards:
    * dam_experiment (title of the pA-DamID experiment, defined in the config: Lamin, H3K9me3,.. )
    * bs (the binsize in the config, now in Kb: 2, 5, 10, ...)
    ---
    Used in the pa_DamID_normalize_track rule.
    """
    bs = wildcards.binsize
    return "{outdir}/pADamID_normalized_tracks/{dam_experiment}-{bs}kb.bw".format(
        bs=int(int(bs) / 1000), **wildcards
    )


###############################################################################
### Rules #####################################################################
###############################################################################


# Rule to generate all the output files
rule all:
    input:
        # all pa-DamID tracks
        expand(
            "{output}/pADamID_results/{site}-{binsize}_{dam_experiment}.txt",
            dam_experiment=config["data_sources"]["pADamID"].keys(),
            site=config["bins"].keys(),
            binsize=[item for sublist in config["bins"].values() for item in sublist],
            output=config["output_directory"],
        ),
        # all coverage files
        expand(
            os.path.join(
                config["output_directory"], "coverage_results/{site}-{bin}_summary.txt"
            ),
            site=config["bins"].keys(),
            bin=[item for sublist in config["bins"].values() for item in sublist],
        ),
        # overlap feature files
        expand(
            os.path.join(
                config["output_directory"],
                "feature_results/overlap_{site}_{feature}.txt",
            ),
            site=config["bins"].keys(),
            feature=config["data_sources"]["feature"].keys(),
            output=config["output_directory"],
        ),
        # nearest feature files
        expand(
            os.path.join(
                config["output_directory"],
                "feature_results/nearest_{site}_{feature}.txt",
            ),
            site=config["bins"].keys(),
            feature=config["data_sources"]["feature"].keys(),
            output=config["output_directory"],
        ),


##################################
### Parsing input ################
##################################


# 1) Convert the summary file to a BED file: one BED per window around the site
# ------------------------------------------------------------------------------
rule expand_integration_sites_to_windows:
    input:
        lambda wildcards: config["summary_file"][wildcards.site],
    output:
        unsorted_bed=temp("{outdir}/inputs/trimmed_{site}-{binsize}.tmp"),
        bed="{outdir}/inputs/trimmed_{site}-{binsize}.bed",
    params:
        genome="mm10",
        size="{binsize}",
        scripts_dir=scripts_dir,
    shell:
        """
        Rscript {scripts_dir}/expand_integrations_to_windows.R \
        --input_bed {input} --output {output.unsorted_bed} --genome {params.genome} \
        --bin_size {params.size}

        bedtools sort -i {output.unsorted_bed} > {output.bed}
        """


##################################
### Coverage information #########
##################################


# 1) Summarize the coverage tracks
# ---------------------------------
# This rule summarizes the coverage tracks for a given site and bin size.
# the output is a single summary file where the columns are the different
# coverage tracks and the rows are the different windows around the site.
# The values are the mean coverage in that window.
rule coverage_summarize_results:
    input:
        bed=lambda wildcards: get_beds(config, wildcards),
        bw=get_coverage_file(config),
    output:
        npz=temp("{outdir}/coverage_results/{site}-{binsize}_summary.npz"),
        txt="{outdir}/coverage_results/{site}-{binsize}_summary.txt",
    params:
        labels=" ".join(get_coverage_ID(config)),
    threads: 5
    shell:
        "multiBigwigSummary BED-file --BED {input.bed}"
        "                            -b {input.bw}"
        "                            -p {threads}"
        "                            --labels {params.labels}"
        "                            --outRawCounts {output.txt}"
        "                            -out {output.npz}"


##################################
### pA-DamID information #########
##################################


# 1) Build the tracks for the pA-DamID data
# ------------------------------------------
# This converts the count files of the pA-DamID into
# the bigwig files for the pA-DamID data and a temporary bedgraph file.
rule paDamID_build_bigwig_track:
    input:
        counts=lambda wildcards: get_paDamID_files(config, wildcards),
        sizes=config["chrom_sizes"]["mm10"],
    output:
        bg=temp("{output}/pADamID_tracks/{dam_experiment}_{condition}_r{rep}.bg"),
        bw="{output}/pADamID_tracks/{dam_experiment}_{condition}_r{rep}-gatc.bw",
    shell:
        "awk -vOFS='\t' '{{"
        "    if (FNR==NR) {{"
        "       cs[$1]=$2"
        "    }} else {{"
        "        mid = ($2+$3)/2;" "        imid = int(mid);"
        "        if (mid!=imid){{"
        "            imid = imid + 1"
        "        }}"
        "        if (FNR==1){{"
        '            printf "%s\t0\t", $1'
        "        }} else if ($2==0) {{"
        '            printf "%.0f\t0\\n%s\t0\t", end, $1'
        "        }}"
        '        printf "%.0f\t0\\n%s\t%.0f\t%.0f\t%.0f\\n%s\t%.0f\t",'
        "               imid-1, $1, imid - 1, imid, $4, $1, imid;"
        "        end=cs[$1]"
        "    }}"
        "}} END {{"
        '    printf "%.0f\t0\\n", end'
        "}}' {input.sizes} <(zcat {input.counts}) | "
        "    sort -k1,1 -k2,2n > {output.bg}; "
        "bedGraphToBigWig {output.bg} {input.sizes} {output.bw}"


# 2) Bin the pA-DamID track
# ------------------------------------------
# From the created bigwig files, we now bin the data into the desired bin size,
# as defined in the config file. The output is a summary file
# (same as the coverage_summarize_results rule), but here every column is a
# different condition of the experiment: experiment or control.
rule paDamID_bin_track:
    input:
        "{output}/pADamID_tracks/{dam_experiment}_experiment_r1-gatc.bw",
        "{output}/pADamID_tracks/{dam_experiment}_experiment_r2-gatc.bw",
        "{output}/pADamID_tracks/{dam_experiment}_control_r1-gatc.bw",
        "{output}/pADamID_tracks/{dam_experiment}_control_r2-gatc.bw",
    output:
        temp("{output}/pADamID_bins/{dam_experiment}-{bs}kb.npz"),
        "{output}/pADamID_bins/{dam_experiment}-{bs}kb.txt",
    params:
        bs=lambda wildcards: int(float(wildcards.bs) * 1000),
        label=["exp1", "exp2", "ctrl1", "ctrl2"],
    threads: 10
    shell:
        "multiBigwigSummary bins --bwfiles {input}"
        "                        --label {params.label}"
        "                        -p {threads}"
        "                        -bs {params.bs}"
        "                        -out {output[0]}"
        "                        --outRawCounts {output[1]}"


# 3) Normalize the pA-DamID track
# ------------------------------------------
# Normalize the binned track by the total number of reads in the track.
# (more datails in the normalize_track.R script)
rule paDamID_normalize_track:
    input:
        counts="{output}/pADamID_bins/{dam_experiment}-{bs}kb.txt",
        sizes=config["chrom_sizes"]["mm10"],
    output:
        "{output}/pADamID_normalized_tracks/{dam_experiment}-{bs}kb.bw",
    params:
        bs=lambda wildcards: int(float(wildcards.bs) * 1000),
        scripts_dir=scripts_dir,
    threads: 5
    shell:
        "Rscript {params.scripts_dir}/normalize_track.R --input {input.counts}"
        "                         --chrom_sizes {input.sizes}"
        "                         --out {output}"
        "                         --binsize {params.bs}"
        "                         --threads {threads}"


# 4) Get the mean pA-DamID signal
# ------------------------------------------
# This rule calculates the mean pA-DamID signal in the given window around the site.
# The output is a text file with the mean signal (both normal and geometric mean)
# for each integration of the site.
#
# In the output file, 'mean' and 'gm_mean' refer to the log2-transformed ratio
# between the experimental data (exp) and the control data (ctrl)
# 'mean_ctrl' and 'gm_mean_ctrl' refer to the log2-transformed control data (ctrl)
rule paDamID_get_means:
    input:
        bed=lambda wildcards: get_beds(config, wildcards),
        counts=lambda wildcards: get_paDamID_counts(config, wildcards),
        bw=lambda wildcards: get_paDamID_normalized_tracks(config, wildcards),
    output:
        "{outdir}/pADamID_results/{site}-{binsize}_{dam_experiment}.txt",
    params:
        window=lambda wildcards: wildcards.binsize,
        scripts_dir=scripts_dir,
    threads: 5
    run:
        input_counts = ",".join(input.counts)
        shell(
            "{params.scripts_dir}/compute_dam_region.R --window {params.window} "
            "                            --bed {input.bed}"
            "                            --bw {input.bw}"
            "                            --counts {input_counts}"
            "                            --out {output}"
            "                            --threads {threads}"
        )


##################################
### nearest feature information ##
##################################
# 0) Expand summary to the specific bin
# ---------------------------------------------------
# This does the same as the first rule here, but now for the nearest feature.
rule expand_integration_sites_to_overlap:
    input:
        lambda wildcards: config["summary_file"][wildcards.site],
    output:
        unsorted_bed=temp("{outdir}/feature_input/input_overlap_{site}.tmp"),
        bed="{outdir}/feature_input/input_overlap_{site}.bed",
    params:
        genome="mm10",
        size=config["bin_for_overlap_feature"],
        scripts_dir=scripts_dir,
    shell:
        """
        Rscript {scripts_dir}/expand_integrations_to_windows.R \
        --input_bed {input} --output {output.unsorted_bed} --genome {params.genome} \
        --bin_size {params.size}

        bedtools sort -i {output.unsorted_bed} > {output.bed}
        """

rule overlap_feature:
    input:
        bed="{outdir}/feature_input/input_overlap_{site}.bed",
        feature=lambda wildcards: config["data_sources"]["feature"][wildcards.feature],
    output:
        "{outdir}/feature_results/overlap_{site}_{feature}.txt",
    shell:
        """
        bedtools intersect -wo -a {input.bed} -b {input.feature} > {output}
        """



rule expand_integration_sites_to_nearest:
    input:
        lambda wildcards: config["summary_file"][wildcards.site],
    output:
        unsorted_bed=temp("{outdir}/feature_input/input_nearest_{site}.tmp"),
        bed="{outdir}/feature_input/input_nearest_{site}.bed",
    params:
        genome="mm10",
        size=config["bin_for_nearest_feature"],
        scripts_dir=scripts_dir,
    shell:
        """
        Rscript {scripts_dir}/expand_integrations_to_windows.R \
        --input_bed {input} --output {output.unsorted_bed} --genome {params.genome} \
        --bin_size {params.size}

        bedtools sort -i {output.unsorted_bed} > {output.bed}
        """

# 1) Get the nearest feature
# --------------------------
# Simply run bedtools closest to get the nearest feature to the integration sites.
# Note that here the summary_file is not expanded to the bins, we only look at the
# exact position the integration happened.
rule nearest_feature:
    input:
        bed="{outdir}/feature_input/input_nearest_{site}.bed",
        feature=lambda wildcards: config["data_sources"]["feature"][wildcards.feature],
    output:
        "{outdir}/feature_results/nearest_{site}_{feature}.txt",
    shell:
        """
        bedtools sort -i {input.feature} | bedtools closest -d -a {input.bed} -b - > {output}
        """