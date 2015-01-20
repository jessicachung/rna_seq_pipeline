#!/bin/env python

'''
Differential expression analysis for single or paired-end RNA-Seq data using 
the Tuxedo protocol (Tophat, Cufflinks) and EdgeR/Voom with counts from 
HTSeq-count.

Authors: Jessica Chung, Bernie Pope, Clare Sloggett, Gayle Philip, Marek Cmero.

Description: This program implements a workflow pipeline for differential
expression analysis. 

It uses the Ruffus library to make the description of the pipeline
more declarative.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

The pipeline is configured by an options file in a python file,
including the actual commands which are run at each stage.
'''

import re
import sys
import copy
from ruffus import *
import os.path
import os
import shutil
from glob import glob
#from rubra.utils import *
from pipeline_base.utils import (runStage, runStageCheck, splitPath, 
                                 getOptions, initLog, getCommand, mkLogFile,
                                 mkTempFilename, getStageOptions, zeroFile,
                                 mkDir)
from pipeline_base.cmdline_args import get_cmdline_args


args = get_cmdline_args()
options = getOptions(args)
logger = initLog(options)


analysis_name = re.sub("\s", "_", options.analysis_name) # Remove any whitespace
samples_csv = options.samples_csv
comparisons_csv = options.comparisons_csv
platform = options.platform
input_dir = options.raw_seq_dir
output_dir = options.output_dir
paired_end = options.paired_end
stranded = options.stranded
trim_reads = options.trim_reads

genome_ref = options.genome_ref
genome_ref_fa = options.genome_ref_fa
gene_ref = options.gene_ref
rrna_ref = options.rrna_ref
cuffdiff_mask_file = options.cuffdiff_mask_file
adapter_seq = options.adapter_seq
seed_mismatches = options.seed_mismatches
palendrome_clip_threshold = options.palendrome_clip_threshold
simple_clip_threshold = options.simple_clip_threshold
trimmomatic_extra_parameters = options.extra_parameters
annotation_dataset = options.annotation_dataset

html_index_script = options.html_index_script
index_script = options.index_script
tophat_script = options.tophat_script
merge_tophat_script = options.merge_tophat_script
fix_tophat_unmapped_reads_script = options.fix_tophat_unmapped_reads_script
htseq_script = options.htseq_script
qc_parse_script = options.qc_parse_script
fastqc_parse_script = options.fastqc_parse_script
alignment_stats_script = options.alignment_stats_script
combine_and_annotate_script = options.combine_and_annotate_script
de_analysis_script = options.de_analysis_script

trimmomatic_path = options.trimmomatic_path
reorder_sam_jar = options.reorder_sam_path
mark_duplicates_jar = options.mark_duplicates_path
rnaseqc_jar = options.rnaseqc_path
add_rg_jar = options.add_or_replace_read_groups_path


java_tmp = "-Djava.io.tmpdir=$TMPDIR" if options.using_merri else ""
samples_csv_name = os.path.basename(samples_csv)[:-4] if \
    samples_csv[-4:].lower() == ".csv" else os.path.basename(samples_csv)
comparisons_csv_name = os.path.basename(comparisons_csv)[:-4] if \
    comparisons_csv[-4:].lower() == ".csv" else \
    os.path.basename(comparisons_csv)


# Output directories:
mkDir(output_dir)
fastqc_dir = os.path.join(output_dir, "fastqc")
mkDir(fastqc_dir)
fastqc_post_trim_dir = os.path.join(output_dir, "fastqc_post_trim")
mkDir(fastqc_post_trim_dir)
transcriptome_dir = os.path.join(output_dir, "transcriptome_index")
# mkDir(transcriptome_dir)
trimmed_dir = os.path.join(output_dir, "trimmed_reads")
mkDir(trimmed_dir)
tophat_raw_dir = os.path.join(output_dir, "tophat_raw")
mkDir(tophat_raw_dir)
tophat_dir = os.path.join(output_dir, "tophat")
mkDir(tophat_dir)
cufflinks_dir = os.path.join(output_dir, "cufflinks")
mkDir(cufflinks_dir)
cuffmerge_dir = os.path.join(output_dir, "cuffmerge")
mkDir(cuffmerge_dir)
cuffdiff_dir = os.path.join(output_dir, "cuffdiff")
mkDir(cuffdiff_dir)
htseq_dir = os.path.join(output_dir, "htseq_count")
mkDir(htseq_dir)
counts_dir = os.path.join(output_dir, "read_counts")
mkDir(counts_dir)
merged_dir = os.path.join(output_dir, "tophat_merged")
mkDir(merged_dir)
rnaseqc_dir = os.path.join(output_dir, "rnaseqc")
mkDir(rnaseqc_dir)
alignment_stats_dir = os.path.join(output_dir, "alignment_stats")
mkDir(alignment_stats_dir)
qc_summary_dir = os.path.join(output_dir, "qc_summary")
mkDir(qc_summary_dir)
main_voom_dir = os.path.join(output_dir, "voom_analysis")
mkDir(main_voom_dir)
main_edger_dir = os.path.join(output_dir, "edgeR_analysis")
mkDir(main_edger_dir)
voom_dir = os.path.join(main_voom_dir, analysis_name + "_voom")
edger_dir = os.path.join(main_edger_dir, analysis_name + "_edgeR")
cuffmerge_sub_dir = os.path.join(cuffmerge_dir, analysis_name + "_cuffmerge")



class sample(object):
    """
    Class for sample information
    """
    def __init__(self, name, condition, covariates = None, files = None):
        self.name = name
        self.condition = condition
        self.files = files if files else []
        self.covariates = covariates if covariates else []
        self.sm = name
    def __repr__(self):
        str = "\n".join(["Name:       %s" % self.name,
                         "Condition:  %s" % self.condition,
                         "Files:      %s" % "\n            ".join(self.files),
                         "Covariates: %s" % self.covariates])
        return str
    def print_info(self):
        str = "\n".join(["\t%s" % self.name,
                         "\t\tCondition:  %s" % self.condition,
                         "\t\tCovariates: %s" % self.covariates,
                         "\t\tFiles:      %s" % \
                         "\n\t\t            ".join(self.files)])
        return str
    def get_trimmed_filenames(self, paired=False):
        if paired:
            return map(lambda x: "%s/%s" % (trimmed_dir, re.sub(r'.fastq.gz',
                       r'.trimmed-paired.fastq.gz', os.path.basename(x))),
                       self.files)
        else:
            return map(lambda x: "%s/%s" % (trimmed_dir, re.sub(r'.fastq.gz',
                       r'.trimmed-single.fastq.gz', os.path.basename(x))),
                       self.files)


def print_heading(heading):
    print "#" * 80 + "\n## %s:" % heading.upper()

print_heading(analysis_name)

# Get files in sequence directory
seqfiles = glob(input_dir + "/*.fastq.gz")
seqfiles.sort()
if not seqfiles:
    print "Error: No *fastq.gz files in sequence directory."
    sys.exit(1)
    
# Open and read CSV files
try:
    with open(samples_csv) as samples_file:
        sample_csv_list = samples_file.read().strip().split("\n")
except IOError:
    print "Error: Cannot open %s" % samples_csv
    sys.exit(1)
try:
    with open(comparisons_csv) as comparison_file:
        comparisons_csv_list = comparison_file.read().strip().split("\n")
except IOError:
    print "Error: Cannot open %s" % comparisons_csv
    sys.exit(1)

# get smrp list from fastq.gz files in sequence directory
smrp_dict = {}
for file in seqfiles:
    try:
        smrp = re.search('(.+\/)?(SM_[A-Za-z0-9-.]+_' \
                         'RP_[A-Za-z0-9-.]+)_?.*', file).group(2)
        if smrp not in smrp_dict:
            smrp_dict[smrp] = [file]
        else:
            smrp_dict[smrp].append(file)
    except AttributeError:
        print "Warning: FASTQ file %s is not in the correct name format. " \
              "File will not be included in analysis." % file

# parse samples in samples.csv
sample_csv_dict = {}
sample_dict = {}
sample_list = []
for i in sample_csv_list:
    if len(i) > 1 and i[0] != "#":
        try:
            line = re.search('([^#]+)#?.+', i).group(1) if "#" in i else i
            line = map(lambda x: x.strip(), line.split(","))
            line[0] = re.search('(SM_)?([A-Za-z0-9-.]+)_?.*', line[0]).group(2)
            if re.search("\+|\-", line[1]):
                print "Error: Non-allowed characters (+,-) in CSV file's " \
                      "condition column."
                sys.exit(1)
            sample_csv_dict[line[0]] = [0] + line[1:]
            sample_dict[line[1]] = []
        except IndexError:
            print "Error: CSV file not formatted correctly.\n" \
                  "Columns in the CSV file should be " \
                  "sample,condition,[covariates...]\n" \
                  "Examples are listed in the user manual document."
            sys.exit(1)

sequences = []
for smrp_name in smrp_dict:
    sm_name = re.search('SM_([A-Za-z0-9-.]+)_RP_.*', smrp_name).group(1)
    try:
        sample_csv_dict[sm_name][0] += 1
        condition = sample_csv_dict[sm_name][1]
        covariates = sample_csv_dict[sm_name][2:]
        smrp_files = smrp_dict[smrp_name]
        new_sample = sample(smrp_name, condition, covariates, smrp_files)
        sample_dict[condition].append(new_sample)
        sample_list.append(new_sample)
        sequences += smrp_files
    except KeyError:
        print "Warning: Sample %s is not listed in the sample CSV file " \
              "and will not be used in analysis." % smrp_name

for sm_name in sample_csv_dict:
    if sample_csv_dict[sm_name][0] == 0:
        print "Warning: Samples %s in CSV file don't match any files in the " \
              "FASTQ directory." % sm_name

if not sequences:
    print "Error: No files given for analysis.\n" \
          "Check if your CSV file and your fastq filenames are in the " \
          "correct format.\n" \
          "Refer to the user manual document for more information."
    sys.exit(1)

# check all samples have the same number of covariates
number_of_covariates = len(sample_list[0].covariates)
for smrp in sample_list:
    if len(smrp.covariates) != number_of_covariates:
        print "Error: Samples in CSV file have an unequal number of " \
              "covariates."
        sys.exit(1)

print_heading("samples")
for condition in sample_dict:
    print "-", condition
    for sample in sample_dict[condition]:
        print sample.print_info()


# parse comparisons in comparisons.csv
comparisons_list = []
comparisons_print = []
no_replicates_warning = False
for i in comparisons_csv_list:
    if len(i) > 1 and i[0] != "#":
        try:
            line = re.search('([^#]+)#?.+', i).group(1) if "#" in i else i
            line = map(lambda x: x.strip(), line.split(","))
            c1 = line[0]
            c2 = line[1]
            n1 = len(sample_dict[c1])
            n2 = len(sample_dict[c2])
            comparisons_list.append([c1, c2])
            comparisons_print.append("\t%s (n = %d) vs. %s (n = %d)" % (c1, n1,
                                                                       c2, n2))
            if n1 < 2 or n2 < 2:
                no_replicates_warning = True
        except IndexError:
            print "Error: CSV file not formatted correctly.\n" \
                  "Columns in the CSV file should be " \
                  "condition1,condition2\n" \
                  "Examples are listed in the user manual document."
            sys.exit(1)
        except KeyError as e:
            print "Error: No samples have condition %s in sample CSV file " \
                  "needed for comparison '%s vs. %s'." % (e, c1, c2)
            sys.exit(1)

print_heading("comparisons")
print "\n".join(comparisons_print)
if no_replicates_warning:
    print "Warning: Lacking replicates for some conditions! " \
          "Some analyses in the pipeline will fail."


@files([samples_csv,
        comparisons_csv],
       ["%s/index.html" % output_dir, 
        "%s/makeIndexHtml.Success" % output_dir])
def makeIndexHtml(inputs, outputs):
    """
    Make index HTML file of results.
    """
    output_filename, flagFile = outputs
    abs_sample_csv_file = os.path.abspath(samples_csv)
    abs_comparison_csv_file = os.path.abspath(comparisons_csv)
    abs_output_dir = os.path.abspath(output_dir)
    runStageCheck('makeIndexHtml', flagFile, logger, options,
                  html_index_script, analysis_name, abs_sample_csv_file,
                  abs_comparison_csv_file, abs_output_dir, output_filename)


if paired_end:
    @transform(sequences,
               regex('(.+\/)?(.+?)\_R1.fastq\.gz'),
               add_inputs(r'\1\2_R2.fastq.gz'), 
               [r'%s/\2_R1_fastqc.zip' % fastqc_dir, 
                r'%s/\2_R2_fastqc.zip' % fastqc_dir, 
                r'%s/\2.fastqc.Success' % fastqc_dir])
    def fastQC(inputs, outputs):
        """
        Obtain stats on reads from fastq using fastQC (paired-end)
        """
        paired1, paired2 = inputs
        out1, out2, flagFile = outputs
        runStageCheck('fastQC', flagFile, logger, options, fastqc_dir, 
                      paired1, paired2)
else:  # if single-end reads
    @transform(sequences, 
               regex('(.+\/)?(.+?)\_R1.fastq\.gz'), 
               [r'%s/\2_R1_fastqc.zip' % fastqc_dir, 
                r'%s/\2.fastqc.Success' % fastqc_dir])
    def fastQC(inputs, outputs):
        """
        Obtain stats on reads from fastq using fastQC (single-end)
        """
        paired1 = inputs
        out1, flagFile = outputs
        paired2 = ""
        runStageCheck('fastQC', flagFile, logger, options, fastqc_dir, 
                      paired1, paired2)


if paired_end and trim_reads:
    @transform(sequences, 
               regex('(.+\/)?(.+?)\_R1.fastq\.gz'),
               add_inputs(r'\1\2_R2.fastq.gz'), 
               [r'%s/\2_R1.trimmed-paired.fastq.gz' % trimmed_dir, 
                r'%s/\2_R1.trimmed-unpaired.fastq.gz' % trimmed_dir, 
                r'%s/\2_R2.trimmed-paired.fastq.gz' % trimmed_dir, 
                r'%s/\2_R2.trimmed-unpaired.fastq.gz' % trimmed_dir, 
                r'%s/\2.trimReads.Success' % trimmed_dir],
               [r'\2'])
    def trimReads(inputs, outputs, samp_name):
        """
        Trim adapter sequences from fastq files using Trimmomatic
        (paired-end)
        """
        paired1, paired2 = inputs
        out1, unpaired1, out2, unpaired2, flagFile = outputs
        paired = "PE"
        parameters = "%s:%d:%d:%d %s" % (adapter_seq, seed_mismatches,
                                         palendrome_clip_threshold,
                                         simple_clip_threshold,
                                         trimmomatic_extra_parameters)
        trim_log = "-trimlog %s/%s.trimReads.log" % \
                (trimmed_dir, samp_name[0]) if options.write_trimmomatic_log \
                else ""
        trimmomatic_input = "%s %s %s %s %s %s ILLUMINACLIP:%s" % \
            (paired1, paired2, out1, unpaired1, out2, unpaired2, parameters)
        runStageCheck('trimReads', flagFile, logger, options, java_tmp,
                      trimmomatic_path, paired, trim_log, trimmomatic_input)
elif trim_reads:
    @transform(sequences, 
               regex('(.+\/)?(.+?)\_R1.fastq\.gz'),
               [r'%s/\2_R1.trimmed-single.fastq.gz' % trimmed_dir, 
                r'%s/\2.trimReads.Success' % trimmed_dir],
               [r'\2'])
    def trimReads(inputs, outputs, samp_name):
        """
        Trim adapter sequences from fastq files using Trimmomatic
        (single-end)
        """
        paired1 = inputs
        out1, flagFile = outputs
        paired = "SE"
        parameters = "%s:%d:%d:%d %s" % (adapter_seq, seed_mismatches,
                                         palendrome_clip_threshold,
                                         simple_clip_threshold,
                                         trimmomatic_extra_parameters)
        trim_log = "-trimlog %s/%s.trimReads.log" % \
                (trimmed_dir, samp_name[0]) if options.write_trimmomatic_log \
                else ""
        trimmomatic_input = "%s %s ILLUMINACLIP:%s" % (paired1, out1,
                                                       parameters)
        runStageCheck('trimReads', flagFile, logger, options, java_tmp, 
                      trimmomatic_path, paired, trim_log, trimmomatic_input)



if paired_end and trim_reads:
    @transform(trimReads, 
               regex('(.+\/)?(.+?)\_R1.trimmed-paired.fastq\.gz'),
               [r'%s/\2_R1.trimmed-paired_fastqc.zip' % fastqc_post_trim_dir,
                r'%s/\2_R2.trimmed-paired_fastqc.zip' % fastqc_post_trim_dir,
                r'%s/\2.fastqcPostTrim.Success' % fastqc_post_trim_dir])
    def fastQCPostTrim(inputs, outputs):
        """
        Obtain stats on reads from fastq using fastQC on trimmed files
        (paired-end)
        """
        paired1, unpaired1, paired2, unpaired2, _success = inputs
        out1, out2, flagFile = outputs    
        runStageCheck('fastQC', flagFile, logger, options,
                      fastqc_post_trim_dir, paired1, paired2)
elif trim_reads:
    @transform(trimReads, 
               regex('(.+\/)?(.+?)\_R1.trimmed-single.fastq\.gz'),
               [r'%s/\2_R1.trimmed-single_fastqc.zip' % fastqc_post_trim_dir,
                r'%s/\2.fastqcPostTrim.Success' % fastqc_post_trim_dir])
    def fastQCPostTrim(inputs, outputs):
        """
        Obtain stats on reads from fastq using fastQC on trimmed files
        (single-end)
        """
        paired1, _success = inputs
        out1, flagFile = outputs
        paired2 = ""
        runStageCheck('fastQC', flagFile, logger, options,
                      fastqc_post_trim_dir, paired1, paired2)

if trim_reads:
    @follows(fastQC)
    @merge(fastQCPostTrim, 
           [r'%s/FastQC_summary.html' % qc_summary_dir,
            r'%s/FastQC_basic_statistics_summary.html' % qc_summary_dir,
            r'%s/FastQC_summary.txt' % qc_summary_dir, 
            r'%s/fastqcSummary.Success' % qc_summary_dir])
    def fastQCSummary(inputs, outputs):
        """
        Parse results from fastQC analysis
        """
        qc_summary, basic_statistics_summary, summary_txt, flagFile = outputs
        paired = "paired" if paired_end else "single"
        runStageCheck('fastQCSummary', flagFile, logger, options, 
                      fastqc_parse_script, fastqc_dir, fastqc_post_trim_dir,
                      qc_summary, basic_statistics_summary, paired, summary_txt)
else:
    @merge(fastQC, 
           [r'%s/FastQC_summary.html' % qc_summary_dir,
            r'%s/FastQC_basic_statistics_summary.html' % qc_summary_dir,
            r'%s/FastQC_summary.txt' % qc_summary_dir, 
            r'%s/fastqcSummary.Success' % qc_summary_dir])
    def fastQCSummary(inputs, outputs):
        """
        Parse results from fastQC analysis
        """
        qc_summary, basic_statistics_summary, summary_txt, flagFile = outputs
        paired = "paired" if paired_end else "single"
        runStageCheck('fastQCSummary', flagFile, logger, options, 
                      fastqc_parse_script, fastqc_dir, fastqc_post_trim_dir,
                      qc_summary, basic_statistics_summary, paired, summary_txt)


@files([genome_ref_fa, gene_ref],
       ["%s/known.rev.1.bt2" % transcriptome_dir, 
        "%s/buildIndex.Success" % transcriptome_dir])
def buildTranscriptomeIndex(inputs, outputs):
    """
    Build index for Bowtie2 from reference files
    """
    genomeRef, geneRef = inputs
    index, flagFile = outputs
    tmp_dir = "%s/tmp" % output_dir
    if not os.path.exists(tmp_dir):
        mkDir(tmp_dir)
    transcriptome_index = "%s/known" % transcriptome_dir
    seq = sequences[0]
    runStageCheck('buildTranscriptomeIndex', flagFile, logger, options, 
                  index_script, seq, tmp_dir, transcriptome_dir, 
                  transcriptome_index, geneRef, genome_ref)



# Get inputs for tophatAlign: 
# Treat samples which have the same SM and RP identifier as technical
#   replicates. Technical replicates are inputted together. 
tophat_files = []
rg_tags = {}
for samp in sample_list:
    output_files = ["%s/%s/accepted_hits.bam" % (tophat_raw_dir, samp.name),
                    "%s/%s.accepted_hits.bam" % (tophat_dir, samp.name),
                    "%s/%s.tophat.Success" % (tophat_raw_dir, samp.name)]
    paired1 = []
    paired2 = []
    if trim_reads:
        trimmed_files = samp.get_trimmed_filenames(paired_end)
        input_files = ["%s/known.rev.1.bt2" % transcriptome_dir] + trimmed_files
        for i in trimmed_files:
            if paired_end:
                if "_R1.trimmed-paired.fastq.gz" in os.path.basename(i):
                    paired1.append(i)
                elif "_R2.trimmed-paired.fastq.gz" in os.path.basename(i):
                    paired2.append(i)
            else:
                if "_R1.trimmed-single.fastq.gz" in os.path.basename(i):
                    paired1.append(i)
                    paired2 = ["False"]
        match = re.search('(.+\/)?SM_([A-Za-z0-9-.]+)_RP_([A-Za-z0-9-.]+)_' \
                          'LB_([A-Za-z0-9-.]+)_ID_([A-Za-z0-9-.]+)_L([0-9]+)_' \
                          'R.\.trimmed-(paired|single)\.fastq\.gz', 
                          input_files[1])
    else:
        untrimmed_files = samp.files
        input_files = ["%s/known.rev.1.bt2" % transcriptome_dir] \
                          + untrimmed_files
        for i in untrimmed_files:
            if paired_end:
                if "_R1.fastq.gz" in os.path.basename(i):
                    paired1.append(i)
                elif "_R2.fastq.gz" in os.path.basename(i):
                    paired2.append(i)
            else:
                if "_R1.fastq.gz" in os.path.basename(i):
                    paired1.append(i)
                    paired2 = ["False"]
        match = re.search('(.+\/)?SM_([A-Za-z0-9-.]+)_RP_([A-Za-z0-9-.]+)_' \
                          'LB_([A-Za-z0-9-.]+)_ID_([A-Za-z0-9-.]+)_L([0-9]+)_' \
                          'R.\.fastq\.gz', input_files[1])
    rgsm = match.group(2) + "_RP-" + match.group(3)   # sample name + replicate
    rglb = match.group(4)                             # library name
    rgid = match.group(5) + "_L" + match.group(6)     # id + lane
    rgpl = platform                                   # platform
    rg_tags[samp.name] = [rgsm, rglb, rgid, rgpl]
    extra_parameters = [",".join(paired1), ",".join(paired2), rgsm, rglb, 
                        rgid, rgpl]
    tophat_files.append([input_files, output_files, extra_parameters])

# print_heading("tophat files")
# for i in tophat_files: print i


if trim_reads:
    @follows(buildTranscriptomeIndex)
    @follows(trimReads)
    @files(tophat_files)
    def tophatAlign(inputs, outputs, extra_parameters):
        """
        Align reads in fastq file using TopHat
        """
        input_files = inputs
        acceptedHits, linkFile, flagFile = outputs
        paired1, paired2, rgsm, rglb, rgid, rgpl = extra_parameters
        sample_dir = os.path.dirname(acceptedHits)
        transcriptome_index = "%s/known" % transcriptome_dir
        runStageCheck('tophatAlign', flagFile, logger, options, tophat_script, 
                      paired1, paired2, sample_dir, gene_ref, genome_ref,
                      transcriptome_index, rgsm, rglb, rgid, rgpl, linkFile)
else:
    @follows(buildTranscriptomeIndex)
    @files(tophat_files)
    def tophatAlign(inputs, outputs, extra_parameters):
        """
        Align reads in fastq file using TopHat
        """
        input_files = inputs
        acceptedHits, linkFile, flagFile = outputs
        paired1, paired2, rgsm, rglb, rgid, rgpl = extra_parameters
        sample_dir = os.path.dirname(acceptedHits)
        transcriptome_index = "%s/known" % transcriptome_dir
        runStageCheck('tophatAlign', flagFile, logger, options, tophat_script, 
                      paired1, paired2, sample_dir, gene_ref, genome_ref,
                      transcriptome_index, rgsm, rglb, rgid, rgpl, linkFile)


@transform(tophatAlign, 
           regex('(.+\/)?(.+?)/accepted_hits\.bam'),
           [r'%s/\2.accepted_hits.sorted.bam' % tophat_dir,
            r'%s/\2.sortBam.Success' % tophat_dir])
def sortBam(inputs, outputs):
    """
    Sort BAM files using Samtools
    """
    originalFile, bamFile, _success = inputs
    output, flagFile = outputs
    output = output[:-4] 
    runStageCheck('sortBam', flagFile, logger, options, bamFile, output)


@transform(sortBam, 
           regex('(.+\/)?(.+?)\.accepted_hits\.sorted\.bam'), 
           [r'%s/\2.accepted_hits.sorted.bam.bai' % tophat_dir, 
            r'%s/\2.indexSortedBam.Success' % tophat_dir])
def indexSortedBam(inputs, outputs):
    """
    Index sorted BAM files using Samtools
    """
    bamFile, _success = inputs
    output, flagFile = outputs
    runStageCheck('indexBam', flagFile, logger, options, bamFile)


@transform(tophatAlign, 
           regex('(.+\/)?(.+?)/accepted_hits\.bam'), 
           add_inputs(r'\1\2/unmapped.bam'), 
           [r'%s/\2.merged.bam' % merged_dir, 
            r'%s/\2.tophatMerge.Success' % merged_dir])
def mergeTophat(inputs, outputs):
    """
    Fix unmapped reads and merges Tophat accepted_hits.bam and unmapped.bam
    """
    [originalFile, bamFile, _success], unmapped = inputs
    output, flagFile = outputs
    sample_dir = os.path.dirname(originalFile)
    runStageCheck('mergeTophat', flagFile, logger, options, 
                  merge_tophat_script, fix_tophat_unmapped_reads_script,
                  sample_dir, output)


@transform(mergeTophat, 
           regex('(.+\/)?(.+?)\.merged\.bam'), 
           [r'%s/\2.merged.reordered.bam' % merged_dir, 
            r'%s/\2.reorderBam.Success' % merged_dir])
def reorderBam(inputs, outputs):
    """
    Reorder BAM files to match the contig ordering of a reference file using
    Picard's Reorder SAM
    """
    bamFile, _success = inputs
    output, flagFile = outputs
    runStageCheck('reorderBam', flagFile, logger, options, java_tmp,
                  reorder_sam_jar, bamFile, output, genome_ref_fa)


@transform(reorderBam, 
           regex('(.+\/)?(.+?).merged\.reordered\.bam'), 
           [r'%s/\2.merged.reordered.addedRG.bam' % merged_dir,
            r'%s/\2.addRG.Success' % merged_dir], [r'\2'])
def addRG(inputs, outputs, samp_name):
    """
    Add Read Groups to BAM file
    """
    bamFile, _success = inputs
    output, flagFile = outputs
    rgsm, rglb, rgid, rgpl = rg_tags[samp_name[0]]
    rgpu = rgid      ### platform unit ?????
    runStageCheck('addRG', flagFile, logger, options, java_tmp, add_rg_jar, 
                  bamFile, output, rgsm, rglb, rgid, rgpl, rgpu)


@transform(addRG, 
           regex('(.+\/)?(.+?)\.merged\.reordered\.addedRG\.bam'), 
           [r'%s/\2.merged.reordered.addedRG.bam.bai' % merged_dir, 
            r'%s/\2.indexReorderedBam.Success' % merged_dir])
def indexReorderedBam(inputs, outputs):
    """
    Index reordered BAM files using Samtools
    """
    bamFile, _success = inputs
    output, flagFile = outputs
    runStageCheck('indexBam', flagFile, logger, options, bamFile)


@follows(indexReorderedBam)
@transform(addRG, 
           regex('(.+\/)?(.+?)\.merged\.reordered\.addedRG\.bam'),
           [r'%s/\2.merged.reordered.addedRG.markdup.bam' % merged_dir,
            r'%s/\2.markdup.log' % merged_dir, 
            r'%s/\2.markDuplicates.Success' % merged_dir])
def markDuplicates(inputs, outputs):
    """
    Mark duplicates in BAM files using Picard
    """
    bamFile, _success = inputs
    output, markDupLog, flagFile = outputs
    runStageCheck('markDuplicates', flagFile, logger, options, java_tmp, 
                  mark_duplicates_jar, bamFile, markDupLog, output)


@transform(markDuplicates, 
           regex('(.+\/)?(.+?)\.merged\.reordered\.addedRG\.markdup\.bam'), 
           [r'%s/\2.merged.reordered.addedRG.markdup.bam.bai' % merged_dir,
            r'%s/\2.indexMardupBam.Success' % merged_dir])
def indexMarkdupBam(inputs, outputs):
    """
    Index marked duplicates BAM files using Samtools
    """
    bamFile, markDupLog, _success = inputs
    output, flagFile = outputs
    runStageCheck('indexBam', flagFile, logger, options, bamFile)


@follows(indexMarkdupBam)
@transform(markDuplicates, 
           regex('(.+\/)?(.+?)\.merged\.reordered\.addedRG\.markdup\.bam'), 
           [r'%s/\2/report.html' % rnaseqc_dir, 
            r'%s/\2.rnaSeQC.Success' % rnaseqc_dir], 
           [r'\2'])
def rnaSeQC(inputs, outputs, samp_name):
    """
    Obtain stats on RNA-seq data using RNA-SeQC
    """
    bamFile, markDupLog, _success = inputs
    output, flagFile = outputs
    samp = "\"%s|%s|%s\"" % (samp_name[0], bamFile, samp_name[0])
    sample_dir = os.path.dirname(output)
    if not rrna_ref:
        rrna = ""
    elif rrna_ref.split(".")[-1] in ("fasta", "fa"):
        rrna = "-BWArRNA %s" % rrna_ref
    elif rrna_ref.split(".")[-1] == "list":
        rrna = "-rRNA %s" % rrna_ref
    else:
        rrna = ""
    paired = "" if paired_end else "-singleEnd"
    runStageCheck('rnaSeQC', flagFile, logger, options, java_tmp, rnaseqc_jar, 
                  paired, samp, genome_ref_fa, gene_ref, rrna, sample_dir)


@follows(indexSortedBam)
@transform(sortBam, 
           regex('(.+\/)?(.+?)\.accepted_hits\.sorted\.bam'),
           [r'%s/\2/transcripts.gtf' % cufflinks_dir,
            r'%s/\2.cufflinksAssembly.Success' % cufflinks_dir])
def cufflinksAssembly(inputs, outputs):
    """
    Assemble aligned reads into transcripts using Cufflinks
    """
    bamFile, _success = inputs
    transcripts, flagFile = outputs    
    samp_dir = os.path.dirname(transcripts)
    runStageCheck('cufflinksAssembly', flagFile, logger, options, samp_dir, 
                  bamFile)


@merge(cufflinksAssembly,
       ['%s/assemblies.txt' % cuffmerge_sub_dir, 
       '%s/assemblies.Success' % cuffmerge_sub_dir])
def createCuffmergeFile(inputs, outputs):
    """
    Create assemblies.txt file containing a list of transcript.gtf files
    """
    assemblies, flagFile = outputs
    transcripts = []
    success = []
    for i in inputs:
        transcripts.append(i[0])
        success.append(i[1])
    mkDir(cuffmerge_sub_dir)
    os.system('echo "%s" > %s' % ("\n".join(transcripts), assemblies))
    os.system('> %s' % flagFile)


@files(createCuffmergeFile, 
       ['%s/merged.gtf' % cuffmerge_sub_dir, 
        '%s/cuffmerge.Success' % cuffmerge_sub_dir])
def cuffmerge(inputs, outputs):
    """
    Create a single merged transcriptome annotation from all assemblies in 
    assembly.txt using cuffmerge
    """
    assemblies, _success = inputs
    transcripts, flagFile = outputs
    runStageCheck('cuffmerge', flagFile, logger, options, gene_ref, 
                  genome_ref_fa, cuffmerge_sub_dir, assemblies)


# Input files in the same group for Cuffdiff analysis
cuffdiff_files = []
for comparison in comparisons_list:
    c1_samples = map(lambda x: x.name, sample_dict[comparison[0]])
    c2_samples = map(lambda x: x.name, sample_dict[comparison[1]])
    c1_files = map(lambda x: "%s/%s.accepted_hits.bam" % (tophat_dir, x),
                   c1_samples)
    c2_files = map(lambda x: "%s/%s.accepted_hits.bam" % (tophat_dir, x),
                   c2_samples)
    label = analysis_name + "_" + "_vs_".join(comparison)
    input_files = ["%s/merged.gtf" % cuffmerge_sub_dir] + c1_files + c2_files
    output_files = ["%s/%s/gene_exp.diff" % (cuffdiff_dir, label),
                    "%s/%s.cuffdiff.Success" % (cuffdiff_dir, label)]
    extra_parameters = [",".join(c1_files), ",".join(c2_files), comparison[0],
                        comparison[1]]
    cuffdiff_files.append([input_files, output_files, extra_parameters])

# print_heading("cuffdiff files")       
# for i in cuffdiff_files: print i


@follows(cuffmerge)
@files(cuffdiff_files)
def cuffdiff(inputs, outputs, extras):
    """
    Identify differentially expressed genes in each group using Cuffdiff.
    """
    merged_gtk = inputs[0]
    output_de, flagFile = outputs
    c1_files, c2_files, c1_label, c2_label = extras
    labels = c1_label + "," + c2_label
    outputDir = os.path.dirname(output_de)
    mask = "-M %s" % cuffdiff_mask_file if cuffdiff_mask_file else ""
    runStageCheck('cuffdiff', flagFile, logger, options, mask, outputDir,
                  labels, merged_gtk, c1_files, c2_files)


@transform(tophatAlign, 
           regex('(.+\/)?(.+?)/accepted_hits\.bam'), 
           [r'%s/\2.accepted_hits.sortedByName.bam' % tophat_dir, 
            r'%s/\2.sortBamByName.Success' % tophat_dir])
def sortBamByName(inputs, outputs):
    """
    Sort BAM file by name
    """
    originalFile, bamFile, _success = inputs
    output, flagFile = outputs
    output = output[:-4] 
    runStageCheck('sortBamByName', flagFile, logger, options, bamFile, output)


@transform(sortBamByName, 
           regex('(.+\/)?(.+?)\.accepted_hits\.sortedByName\.bam'), 
           add_inputs(r'%s/\2/unmapped.bam' % tophat_raw_dir), 
           [r'%s/\2.alignmentStats.txt' % alignment_stats_dir, 
            r'%s/\2.alignmentStats.Success' % alignment_stats_dir])
def alignmentStats(inputs, outputs):
    """
    Count the number of reads which had unique alignments, the number of
    reads which had multiple alignments, and the number of unmapped reads
    """
    [bamFile, _success], unmappedBam = inputs
    output, flagFile = outputs
    paired = "paired" if paired_end else "single"
    runStageCheck('alignmentStats', flagFile, logger, options, 
                  alignment_stats_script, bamFile, unmappedBam, output, paired)


@follows(alignmentStats)
@follows(fastQCSummary)
@merge(rnaSeQC, 
       [r'%s/qc_summary.html' % qc_summary_dir, 
        r'%s/qcSummary.Success' % qc_summary_dir])
def qcSummary(inputs, outputs):
    """
    Parse results from QC analysis
    """
    qc_summary, flagFile = outputs
    paired = "paired" if paired_end else "single"
    runStageCheck('qcSummary', flagFile, logger, options, qc_parse_script, 
                  fastqc_dir, fastqc_post_trim_dir, alignment_stats_dir,
                  rnaseqc_dir, qc_summary, paired)


@transform(sortBamByName, 
           regex('(.+\/)?(.+?)\.accepted_hits\.sortedByName\.bam'),
           [r'%s/\2.union_HTSeqCount.txt' % htseq_dir,
            r'%s/\2.strictIntersect_HTSeqCount.txt' % htseq_dir,
            r'%s/\2.htseqCount.Success' % htseq_dir])
def countReads(inputs, outputs):
    """
    Count reads for each feature in GTF file.
    """
    bamFile, _success = inputs
    unionFile, strictFile, flagFile = outputs
    runStageCheck('countReads', flagFile, logger, options, htseq_script, 
                  bamFile, gene_ref, unionFile, strictFile, stranded)


@merge(countReads, 
       ['%s/%s_samples.csv' % (counts_dir, analysis_name),
        '%s/%s_comparisons.csv' % (counts_dir, analysis_name),
        r'%s/%s_counts.txt' % (counts_dir, analysis_name),
        r'%s/%s_counts.RData' % (counts_dir, analysis_name),
        r'%s/%s_counts.stdout' % (counts_dir, analysis_name),
        r'%s/%s_counts.stderr' % (counts_dir, analysis_name),
        r'%s/%s.combineAndAnnotate.Success' % (counts_dir, analysis_name)])
def combineAndAnnotate(inputs, outputs):
    """
    Create csv files containing sample information and comparison information
    needed for edgeR and voom analysis. Combine feature counts from HTSeq into
    one CSV file. Also removes all features with zero counts. Annotates Ensembl
    IDs with symbols, chr, description etc. Text file of raw counts can be used
    for DGE-Vis.
    """
    sample_R_csv, comparison_R_csv, plain_text_counts, rdata_counts, \
            combine_stdout, combine_stderr, flagFile = outputs
    # If replicates labels are all identical, then remove RP tag
    rp_list = map(lambda x: re.search('SM_[A-Za-z0-9-.]+_RP_([A-Za-z0-9-.]+)',
                                      x.name).group(1), sample_list)
    if len(set(rp_list)) == 1:
        for s in sample_list:
            s.sm = re.search('SM_([A-Za-z0-9-.]+)_RP_([A-Za-z0-9-.]+)',
                             s.name).group(1)
    # Write sample csv file
    try:
        with open(sample_R_csv, 'w') as output_file:
            output_lines = []
            for smrp_name in sample_list:
                htseq_count = "%s/%s.union_HTSeqCount.txt" % (htseq_dir, 
                                                              smrp_name.name)
                if smrp_name.covariates:
                    output_lines.append(",".join([smrp_name.sm, htseq_count,
                            smrp_name.condition, 
                            ",".join(smrp_name.covariates)]))
                else:
                    output_lines.append(",".join([smrp_name.sm, htseq_count,
                            smrp_name.condition]))
            output_file.write("\n".join(output_lines) + "\n")
    except:
        print "Error. Could not create file %s" % sample_R_csv
        sys.exit(1)
    # Write comparison csv file
    try:
        with open(comparison_R_csv, 'w') as output_file:
            output_lines = []
            for comparison in comparisons_list:
                output_lines.append(",".join(comparison))
            output_file.write("\n".join(output_lines) + "\n")
    except:
        print "Error. Could not create file %s" % comparison_R_csv
        sys.exit(1)
    annotation_dataset = str(options.annotation_dataset)
    runStageCheck('combineAndAnnotate', flagFile, logger, options, 
                  sample_R_csv, comparison_R_csv, plain_text_counts,
                  rdata_counts, annotation_dataset, 
                  combine_and_annotate_script, combine_stdout, combine_stderr)


@files(combineAndAnnotate, 
       ['%s/voom.stdout' % voom_dir, 
        '%s/voom.stderr' % voom_dir, 
        '%s/voom.Success' % voom_dir])
def voom(inputs, outputs):
    """
    Perform DE analysis using Voom (limma)
    """
    sample_R_csv, comparison_R_csv, plain_text_counts, rdata_counts, \
            combine_stdout, combine_stderr, _success = inputs
    voom_stdout, voom_stderr, flagFile = outputs
    mkDir(voom_dir)
    runStageCheck('voom', flagFile, logger, options, rdata_counts, voom_dir,
                  "voom", de_analysis_script, voom_stdout, voom_stderr)


@files(combineAndAnnotate, 
       ['%s/edgeR.stdout' % edger_dir, 
        '%s/edgeR.stderr' % edger_dir, 
        '%s/edgeR.Success' % edger_dir])
def edgeR(inputs, outputs):
    """
    Perform DE analysis using edgeR
    """
    sample_R_csv, comparison_R_csv, plain_text_counts, rdata_counts, \
            combine_stdout, combine_stderr, _success = inputs
    edger_stdout, edger_stderr, flagFile = outputs
    mkDir(edger_dir)
    runStageCheck('edgeR', flagFile, logger, options, rdata_counts, edger_dir,
                  "edgeR", de_analysis_script, edger_stdout, edger_stderr)




# Invoke the pipeline.
pipelineOptions = options.pipeline
endTasks = pipelineOptions['end']
forcedTasks = pipelineOptions['force']
style = pipelineOptions['style']

if style == 'run':
    # Perform the pipeline steps.
    pipeline_run(endTasks, 
                 multiprocess = pipelineOptions['procs'],
                 logger = black_hole_logger, 
                 forcedtorun_tasks = forcedTasks,
                 gnu_make_maximal_rebuild_mode = options.maximal_rebuild_mode)
elif style == 'flowchart':
    # Draw the pipeline as a diagram.
    pipeline_printout_graph('flowchart.svg', 'svg', endTasks,
                            no_key_legend = False)
elif style == 'print':
    pipeline_printout(sys.stdout, endTasks, verbose = 5, wrap_width=100000, 
            forcedtorun_tasks = forcedTasks,
            gnu_make_maximal_rebuild_mode = options.maximal_rebuild_mode)



