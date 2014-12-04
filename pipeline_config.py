#---------------------------------
# PIPELINE RUN
#---------------------------------
# The configuration settings to run the pipeline. These options are overwritten
# if a new setting is specified as an argument when running the pipeline.
# These settings include:
# - logDir: The directory where the batch queue scripts are stored, along with 
#   stdout and stderr dumps after the job is run.
# - logFile: Log file in logDir which all commands submitted are stored.
# - style: the style which the pipeline runs in. One of:
#    - 'print': prints the stages which will be run to stdout,
#    - 'run': runs the pipeline until the specified stages are finished, and 
#    - 'flowchart': outputs a flowchart of the pipeline stages specified and
#       their dependencies.
#  - procs: the number of python processes to run simultaneously. This 
#    determines the  maximum parallelism of the pipeline. For distributed jobs
#    it also constrains the maximum total jobs submitted to the queue at any one
#    time.
#  - verbosity: one of 0 (quiet), 1 (normal), 2 (chatty).
#  - end: the desired tasks to be run. Rubra will also run all tasks which are
#    dependencies of these tasks.
#  - force: tasks which will be forced to run, regardless of timestamps.
#  - rebuild: one of 'fromstart','fromend'. Whether to calculate which 
#    dependencies will be rerun by working back from an end task to the latest
#    up-to-date task, or forward from the earliest out-of-date task. 'fromstart'
#    is the most conservative and commonly used as it brings all intermediate 
#    tasks up to date.
#  - manager: "pbs" or "slurm"

pipeline = {
    "logDir": "log",
    "logFile": "pipeline_commands.log",
    "style": "print",
    "procs": 16,
    "verbose": 2,
    "end": ["fastQCSummary", "voom", "edgeR", "qcSummary"],
    "force": [],
    "rebuild": "fromstart",
    "manager": "slurm",
}


# This option specifies whether or not you are using VLSCI's Merri or Barcoo
# cluster. If True, this changes java's tmpdir to the job's tmp dir on 
# /scratch ($TMPDIR) instead of using the default /tmp which has limited space.
using_merri = True


# Optional parameter governing how Ruffus determines which part of the 
# pipeline is out-of-date and needs to be re-run. If set to False, Ruffus
# will work back from the end target tasks and only execute the pipeline
# after the first up-to-date tasks that it encounters.
# Warning: Use with caution! If you don't understand what this option does, 
# keep this option as True.
maximal_rebuild_mode = True


#---------------------------------
# CONFIG
#---------------------------------

# Name of analysis. Changing the name will create new sub-directories for
# voom, edgeR, and cuffdiff analysis.
analysis_name = "analysis_v1"

# The directory containing *.fastq.gz read files.
raw_seq_dir = "/path_to_project/fastq_files/"

# Path to the CSV file with sample information regarding condition and
# covariates if available.
samples_csv = "/path_to_project/fastq_files/samples.csv"

# Path to the CSV file with which comparisons to make.
comparisons_csv = "/path_to_project/fastq_files/comparisons.csv"

# The output directory.
output_dir = "/path_to_project/results/"

# Sequencing platform for read group information.
platform = "Illumina"

# If the experiment is paired-end or single-end: True (PE) or False (SE).
paired_end = False

# Whether the experiment is strand specific: "yes", "no", or "reverse".
stranded = "no"


#---------------------------------
# REFERENCE FILES
#---------------------------------
# Most reference files can be obtained from the Illumina iGenomes project:
# http://cufflinks.cbcb.umd.edu/igenomes.html

# Bowtie 2 index files: *.1.bt2, *.2.bt2, *.3.bt2, *.4.bt2, *.rev.1.bt2, 
# *.rev.2.bt2.
genome_ref = "/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bowtie_Indexed/human_g1k_v37"

# Genome reference FASTA. Also needs an indexed genome (.fai) and dictionary 
# (.dict) file in the same directory.
genome_ref_fa = "/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bowtie_Indexed/human_g1k_v37.fa"

# Gene set reference file (.gtf). Recommend using the GTF file obtained from 
# Ensembl as Ensembl gene IDs are used for annotation (if specified).
gene_ref = "/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/TuxedoSuite_Ref_Files/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"

# Either a rRNA reference fasta (ending in .fasta or .fa) or an GATK interval 
# file (ending in .list) containing rRNA intervals to calculate the rRNA 
# content. Can set as False if not available.
# rrna_ref = "/vlsci/VR0002/shared/Reference_Files/rRNA/human_all_rRNA.fasta"
rrna_ref = "/vlsci/VR0002/shared/jchung/human_reference_files/human_rRNA.list"

# Optional tRNA and rRNA sequences to filter out in Cuffdiff (.gtf or .gff). 
# Set as False if not provided.
cuffdiff_mask_file = False


#---------------------------------
# TRIMMOMATIC PARAMETERS
#---------------------------------
# Parameters for Trimmomatic (a tool for trimming Illumina reads).
# http://www.usadellab.org/cms/index.php?page=trimmomatic

# Path of a FASTA file containing adapter sequences used in sequencing.
adapter_seq = "/vlsci/VR0002/shared/jchung/human_reference_files/TruSeqAdapters.fa"

# The maximum mismatch count which will still allow a full match to be 
# performed.
seed_mismatches = 2

# How accurate the match between the two 'adapter ligated' reads must be for 
# PE palindrome read alignment.
palendrome_clip_threshold = 30

# How accurate the match between any adapter etc. sequence must be against a 
# read.
simple_clip_threshold = 10

# The minimum quality needed to keep a base and the minimum length of reads to 
# be kept.
extra_parameters = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Output Trimmomatic log file
write_trimmomatic_log = True

#---------------------------------
# R PARAMETERS
#---------------------------------

# Get annotations from Ensembl BioMart. GTF file needs to use IDs from Ensembl. 
# Set as False to skip annotation, else 
# provide the name of the dataset that will be queried. Attributes to be
# obtained include gene symbol, chromosome name, description, and gene biotype.
# Commonly used datasets:
# human: "hsapiens_gene_ensembl"
# mouse: "mmusculus_gene_ensembl"
# rat: "rnorvegicus_gene_ensembl"
# You can list all available datasets in R by using the listDatasets fuction:
#  > library(biomaRt)
#  > listDatasets(useMart("ensembl"))
# The gene symbol is obtained from the attribute "hgnc_symbol" (human) or 
# "mgi_symbol" (mice/rats) if available. If not, the "external_gene_id" is used
# to obtain the gene symbol. You can change this by editing the script:
# scripts/combine_and_annotate.r

annotation_dataset = "hsapiens_gene_ensembl"


#---------------------------------
# SCRIPT PATHS
#---------------------------------
# Paths to other wrapper scripts needed to run the pipeline. Make sure these 
# paths are relative to the directory where you plan to run the pipeline in or
# change them to absolute paths.

html_index_script = "scripts/html_index.py"
index_script = "scripts/build_index.sh"
tophat_script = "scripts/run_tophat.sh"
merge_tophat_script = "scripts/merge_tophat.sh"
fix_tophat_unmapped_reads_script = "scripts/fix_tophat_unmapped_reads.py"
htseq_script = "scripts/run_htseq.sh"
fastqc_parse_script = "scripts/fastqc_parse.py"
qc_parse_script = "scripts/qc_parse.py"
alignment_stats_script = "scripts/alignment_stats.sh"
combine_and_annotate_script = "scripts/combine_and_annotate.R"
de_analysis_script = "scripts/de_analysis.R"


#---------------------------------
# PROGRAM PATHS
#---------------------------------
trimmomatic_path = "/usr/local/trimmomatic/0.30/trimmomatic-0.30.jar"
reorder_sam_path = "/usr/local/picard/1.69/lib/ReorderSam.jar"
mark_duplicates_path = "/usr/local/picard/1.69/lib/MarkDuplicates.jar"
rnaseqc_path = "/usr/local/rnaseqc/1.1.7/RNA-SeQC_v1.1.7.jar"
add_or_replace_read_groups_path = "/usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar"

