#!/usr/bin/python3.6
"""
Sequence processing pipeline used to analyze packaged viral genomes
"""
from argparse import ArgumentParser
from csv import DictReader, DictWriter
import logging
import os
from pathlib import Path
from subprocess import run
import sys
import csv

# Argument parsing setup
parser = ArgumentParser(description='Process sequencing files '
                                    'and collect stats')

parser.add_argument('base_dir', type=str,
                    help='Base of where processing takes place. All paths '
                         'in the csv are assumed to be relative to this path '
                         'and results will be placed in "Output" directory '
                         'within this path.')

parser.add_argument('csv_file', type=str,
                    help='CSV file detailing samples, and where the relevant '
                         'files for those samples can be found, all paths '
                         'are relative to the base_dir.')

parser.add_argument('AF', type=float,
					help='Allele fraction: the fraction of all reads conferring '
						 'a base. AF must be above this threshold for pipeline to ' 
						 'confidently call a base at a specified position.')

# Logging setup
log = logging.getLogger("fastseq")
log.addHandler(logging.StreamHandler(sys.stdout))

logfile_formatter = logging.Formatter('%(asctime)s - %(name)s - '
                                      '%(levelname)s - %(message)s')
logfile_handler = logging.FileHandler('fastseq.log')
logfile_handler.setFormatter(logfile_formatter)
log.addHandler(logfile_handler)

log.setLevel(logging.INFO)

args = parser.parse_args()

BASE_DIR = Path(args.base_dir)
CSV_PATH = Path(args.csv_file)
STATS_OUTPUT_PATH = BASE_DIR / "Output" / "final_stats.csv"

# ----------------------------------
# Configuration Variable Definitions
# ----------------------------------

# For simplicity, all program configurations + paths are treated as global
# These are designed to work with the corresponding docker image, can be
# tweaked to work in different contexts.

# Tool Paths
TRIMMOMATIC = "/tools/trimmomatic/trimmomatic-0.38.jar"
BWA = "/tools/bwa/bwa"
SAMTOOLS = "/tools/samtools/bin/samtools"
BCFTOOLS = "/tools/samtools/bin/bcftools"
PICARD = "/tools/picard/picard.jar"

GATK = "/Users/lucy.maynard/Downloads/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"

# Configuration for Trimmomatic
LEAD_SCORE = 3
TRAIL_SCORE = 3
MIN_LEN = 50
WINDOW_SIZE = 4
WINDOW_QUALITY = 20

# Configuration for bcftools
VCF_QUAL = 20
VCF_DP = 10
VCF_AF = AF

# Configuration for Picard
PICARD_COVERAGE_CAP = 100000
PICARD_FAST_ALG = True
PICARD_SAMPLE_SIZE = 5000


# ----------------------------------
# Function Definitions
# ----------------------------------

def configure_paths(sample, fwd_read, rev_read, adapter_pth, ref_pth):
    """
    Create all derived paths based on fwd/rev read, adapter, reference
    Also sets up an output directory with sample name to output files
    Args:
        sample (str): Sample name
        fwd_read (str): Path to forward read rel. to docker base in .fastq.gz
        rev_read (str): Path to forward read rel. to docker base in .fastq.gz
        adapter_pth (str): Path to adapter rel. to docker base in .fasta
            see trimmomatic documentation for what to name the sequences in the
            .fasta file
        ref_pth (str): Path to reference rel. to docker base in .fasta
    Returns:
        dict: A dictionary with keys of type str, values of type Path,
        See function for what keys map to what.
    """

    sample_base = BASE_DIR / "Output" / sample
    ref_dict = ref_pth.split(".fasta")[0] + str(".dict")
    print(ref_dict)
    os.makedirs(sample_base)

    return {
        "output_base": sample_base,

        "fwd_read": BASE_DIR / fwd_read,
        "rev_read": BASE_DIR / rev_read,
        "adapter_pth": BASE_DIR / adapter_pth,
        "ref_pth": BASE_DIR / ref_pth,
        "ref_dict": BASE_DIR / ref_dict,

        # Derived Sample Paths
        "fwd_trimmed": BASE_DIR / f"{fwd_read}.trimmed.fastq",
        "rev_trimmed": BASE_DIR / f"{rev_read}.trimmed.fastq",

        "fwd_unpaired": BASE_DIR / f"{fwd_read}.unpaired.fastq",
        "rev_unpaired": BASE_DIR / f"{rev_read}.unpaired.fastq",

        "sam_file": sample_base / f"{sample}.sam",
        "bam_file": sample_base / f"{sample}.bam",

        "consensus_file": sample_base / f"{sample}.consensus.fasta",
        "readgroup_bam": sample_base / f"{sample}.readgroup.bam",

        "filtered_vcf_file": sample_base / f"{sample}.filtered.vcf.gz",
        "vcf_stats_file": sample_base / f"{sample}.vcf.stats.txt",

        "wgs_metrics_file": sample_base / f"{sample}.picard_wgs.txt",
        "size_metrics_file": sample_base / f"{sample}.picard_size.txt",
        "size_histogram_file": sample_base / f"{sample}.picard_size_hist.pdf"
    }


def trimmomatic(sample, paths):
    """
    Simple wrapper for applying trimmomatic, trims adapters and cleans
    sequence ends. Uses phred33 quality threshold.

    Args:
        sample (str): Name of sample
        paths (dict): Paths collection

    Returns: None
    """

    log.info(f"Starting trimmomatic for {sample}...")

    run(["java", "-jar", TRIMMOMATIC, "PE", "-phred33",
         paths["fwd_read"], paths["rev_read"],  # Input Files

         # Output Files
         paths["fwd_trimmed"], paths["fwd_unpaired"],
         paths["rev_trimmed"], paths["rev_unpaired"],

         f"ILLUMINACLIP:{paths['adapter_pth']}:4:20:10",
         f"LEADING:{LEAD_SCORE}",
         f"TRAILING:{TRAIL_SCORE}",
         f"SLIDINGWINDOW:{WINDOW_SIZE}:{WINDOW_QUALITY}",
         f"MINLEN:{MIN_LEN}"])

    log.info(f"...end trimmomatic for {sample}.")


def bwa(sample, paths):
    """
    Simple wrapper for applying BWA. First indexes then applys mem algorithm

    Args:
        sample (str): Name of sample
        paths (dict): Paths collection

    Returns: None
    """

    # index reference
    log.info(f"Starting BWA Index for {sample}...")

    run([BWA, "index", paths["ref_pth"]])

    log.info(f"...end BWA Index for {sample}.")

    # mem algorithm to align reads + generate .sam file
    log.info(f"Starting BWA mem for {sample}...")

    with open(paths["sam_file"], "w") as f:
        run([BWA, "mem",
             paths["ref_pth"], paths["fwd_trimmed"], paths["rev_trimmed"]],
            stdout=f)  # output to SAMPLE_SAM_PTH for samtools

    log.info(f"...end BWA mem for {sample}.")


def samtools_gatk(sample, paths):
    """
    Wrapper for applying samtools/bcftools/GATK for indel/SNP calling +
    consensus sequence generation.

    First converts SAM file to BAM format, then generates a reference
    dictionary, and repairs readgroups in BAM file, which is required for GATK.
    Finally creates a VCF file and filters it by allele depth & fraction
    Args:
        sample (str): Name of sample
        paths (dict): Paths collection
    Returns: None
    """

    # convert .sam to .bam
    log.info(f"Starting conversion of sam to bam for {sample}...")

    with open(paths["bam_file"], "w") as f:
        run([SAMTOOLS, "sort", paths["sam_file"]],
            stdout=f)  # output to SAMPLE_BAM_PTH
    run([SAMTOOLS, "index", paths["bam_file"]])
    run([SAMTOOLS, "faidx", paths["ref_pth"]])

    log.info(f"...end conversion of sam to bam for {sample}...")

    # generate reference dictionary
    log.info(f"Creating reference dictionary for {sample}...")

    run(["java", "-Xmx2048m", "-jar", GATK, "CreateSequenceDictionary",
         "-R", paths["ref_pth"],
         "-O", paths["ref_dict"]])  # Output file

    log.info(f"...done creating reference dictionary for {sample}...")

    # Fix read groups in bam file
    log.info(f"Fixing read groups for {sample}...")

    run(["java", "-Xmx2048m", "-jar", GATK, "AddOrReplaceReadGroups",
         "-I", paths["bam_file"],  # Input file
         "-O", paths["readgroup_bam"],  # Reference file
         "-RGLB", "lib1",
         "-RGPL", "illumina",
         "-RGPU", "unit1",
         "-RGSM", "20"])  # Output file
    run([SAMTOOLS, "index", paths["readgroup_bam"]])

    log.info(f"...done fixing read groups for {sample}...")

    # generate VCF file with GATK Haplotype Caller
    log.info(f"Creating VCF file for {sample}...")

    run(["java", "-Xmx2048m", "-jar", GATK, "HaplotypeCaller",
         "-R", paths["ref_pth"],  # Input file
         "-I", paths["readgroup_bam"],  # Reference file
         "-O", paths["vcf_file"]])  # Output file

    log.info(f"...end creating VCF file for {sample}...")

    # filter .vcf file by quality thresholds and output to a compressed file
    log.info(f"Starting VCF filter for {sample}...")

    with open(paths["filtered_vcf_file"], "w") as f:
        run([BCFTOOLS, "filter", "-i",
             f"FORMAT/DP>{VCF_DP} && "
             f"(FORMAT/AD[*:1]/ FORMAT/DP)>{VCF_AF} && "
             "FORMAT/AD[*:1] != '*'",
             "-Oz",
             paths["vcf_file"]], stdout=f)

    log.info(f"...end VCF filter for {sample}.")

    # index compressed and filtered VCF file
    run(["tabix", paths["filtered_vcf_file"]])

    # create consensus sequence
    log.info(f"Starting consensus generation for {sample}...")

    with open(paths["consensus_file"], "w") as f:
        run([BCFTOOLS, "consensus", "-f", paths["ref_pth"],
             paths["filtered_vcf_file"]], stdout=f)

    log.info(f"...end consensus filter for {sample}.")


def generate_stats(sample, paths):
    """
    Wrapper to compute stats from bcf tools and from picard.
    Gets VCF stats from BCF tools then collects WGS and Size metrics using
    Picard.
    Args:
        sample (str): Name of sample
        paths (dict): Paths collection
    Returns: None
    """
    # BCFTOOLS stats
    log.info(f"Starting VCF stats for {sample}...")

    with open(paths["vcf_stats_file"], "w") as f:
        run([BCFTOOLS, "stats", paths["filtered_vcf_file"]],
            stdout=f)  # output to SAMPLE_VCF_STATS_PTH

    log.info(f"...end VCF stats for {sample}.")

    # Picard CollectWgsMetrics (library alignment stats)
    log.info(f"Starting picard WGS stats for {sample}...")

    run(["java", "-Xmx2048m", "-jar", PICARD, "CollectWgsMetrics",
         f"COVERAGE_CAP={PICARD_COVERAGE_CAP}",
         f"USE_FAST_ALGORITHM={PICARD_FAST_ALG}",
         f"SAMPLE_SIZE={PICARD_SAMPLE_SIZE}",
         f"I={paths['bam_file']}",  # Input file
         f"R={paths['ref_pth']}",  # Reference file
         f"O={paths['wgs_metrics_file']}"])  # Output file

    log.info(f"...end picard WGS stats for {sample}.")

    # Picard CollectInsertSizeMetrics (fragment size stats)
    log.info(f"Starting picard size stats for {sample}...")

    run(["java", "-Xmx2048m", "-jar", PICARD, "CollectInsertSizeMetrics",
         f"I={paths['bam_file']}",
         f"H={paths['size_histogram_file']}",
         f"O={paths['size_metrics_file']}"])

    log.info(f"...end picard size stats for {sample}.")


def extract_bcf_stats(path):
    """
    Extract relevant information from BCF Stats file for a single sample

    Specifically extract SNPs, MNPs, indels, "others", multiallelic sites,
    and multiallelic SNPsites.

    No effort is made to convert strings to numbers for the stat values.

    Args:
        path (str): path to the BCF stats file

    Returns:
        dict: keys as stat names and the values as stat values
    """

    # Not ideal to hardcode here nor below, but gets the job done
    stats_of_interest = {"number of SNPs:",
                         "number of MNPs:",
                         "number of indels:",
                         "number of others:",
                         "number of multiallelic sites:",
                         "number of multiallelic SNP sites:"}
    p = str(path.resolve()).split("/")[-1]
    p2 = p.split(".")[0]
    stats = {"Sample": p2}

    with open(path) as statsf:
        for line in statsf:

            if line.startswith("SN"):
                parts = line.strip().split("\t")
                stat = parts[-2]
                num = parts[-1]

                if stat in stats_of_interest:
                    stats[stat.strip(":")] = num

    return stats


def extract_picard_stats(path):
    """
    Extract relevant information from picard wgs or size stats file.

    This is assumed to be for a single sample and that there will only be
    two lines in the "METRICS CLASS" section, which is the only section we'll
    extract.

    No effort is made to convert strings to numbers for the stat values.

    Args:
        path (str): path to the picard wgs stats file

    Returns:
        dict: keys as stat names and the values as stat values
    """

    with open(path) as statsf:
        split_lines = []
        keep_line = False
        for line in statsf:
            if keep_line:
                split_lines.append(line.strip().split("\t"))

            # if we see metrics label, set flag to start collecting data
            if line.startswith("## METRICS CLASS"):
                keep_line = True

            # stop at first empty line, though in practice we expect this
            # to happen after exactly 2 lines read
            if keep_line and not line.strip():
                break

    # expecting only 2 lines, header row and values row
    stats = dict(zip(split_lines[0], split_lines[1]))

    return stats


# ----------------------------------
# Main Code Execution
# ----------------------------------

with open(CSV_PATH) as csvfile:
    reader = DictReader(csvfile)

    final_stats = []

    for entry in reader:
        sample_name = entry["Sample"]
        fwd_pth = entry["Forward Read Path"]
        rev_pth = entry["Reverse Read Path"]
        ad_pth = entry["Adapter Path"]
        rf_pth = entry["Reference Path"]

        path_dict = configure_paths(sample_name, fwd_pth, rev_pth,
                                    ad_pth, rf_pth)

        # 1. Trimmomatic (trim adapters and filter by quality threshold) PE
        # (paired end algorithm) with -phred33 (quality threshold)
        trimmomatic(sample_name, path_dict)

        # 2. BWA (align to reference)
        bwa(sample_name, path_dict)

        # 3. SAMTOOLS/BCFTOOLS (call SNPS/indels)
        samtools_gatk(sample_name, path_dict)

        # 4. Generate statistics
        generate_stats(sample_name, path_dict)

        # 5. Extract statistics and collate into a single row
        vcf_st = extract_bcf_stats(path_dict["vcf_stats_file"])
        picard_wgs_st = extract_picard_stats(path_dict["wgs_metrics_file"])
        picard_size_st = extract_picard_stats(path_dict["size_metrics_file"])

        # Assuming no overlap in stat names
        vcf_st.update(picard_wgs_st)
        vcf_st.update(picard_size_st)

        final_stats.append(vcf_st)

log.info(f"Starting writing final stats...")

with open(STATS_OUTPUT_PATH, "w") as statsout:
    # Extract keys (Metric names, e.g. number of snps) from first dictionary
    # of final_stat list
    values_as_cols = [[*final_stats[0]]]
    # Extract values from each dict in final_stats list and append as elemnt
    # to values_as_cols list
    for item in final_stats:
        vals = list(item.values())
        values_as_cols.append(vals)
    values_as_rows = list(zip(*values_as_cols))
    compiled_stats = [list(elem) for elem in values_as_rows]

    # Write compiled_stats list (of lists) to csv file
    writer = csv.writer(statsout, delimiter=",")
    for item in compiled_stats:
        writer.writerow(item)

log.info(f"...end writing stats.")
