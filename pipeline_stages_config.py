#---------------------------------
# BATCH DEFAULTS
#---------------------------------
# The default options applied to each stage in the shell script. These options
# are overwritten by the options provided by the individual stages in the next
# section.
# Stage options which Rubra will recognise are: 
#  - distributed: a boolean determining whether the task should be submitted to
#    a cluster job scheduling system (True) or run on the system local to Rubra
#    (False).   
#  - walltime: for a distributed PBS job, gives the walltime requested from the
#    job queue system; the maximum allowed runtime. For local jobs has no
#    effect.
#  - memInGB: for a distributed PBS job, gives the memory in gigabytes requested
#    from the job queue system. For local jobs has no effect.
#  - queue: for a distributed PBS job, this is the name of the queue to submit
#    the job to. For local jobs has no effect.
#  - modules: the modules to be loaded before running the task. This is intended
#    for systems with environment modules installed. Rubra will call module load
#    on each required module before running the task. Note that defining modules
#    for individual stages will override (not add to) any modules listed here.
#    This currently only works for distributed jobs.

stageDefaults = {
    "distributed": True,
    "walltime": "01:00:00",
    "memInGB": 4,
    "queue": "main", 
    "modules": [
        "perl/5.18.0",
        "java/1.7.0_25",
        "samtools-intel/0.1.19",
        "python-gcc/2.7.5",
        "fastqc/0.10.1",
        "bowtie2-intel/2.1.0",
        "tophat-gcc/2.0.8",
        "cufflinks-gcc/2.1.1",
        "bwa-intel/0.7.5a",
    ],
    "manager": "slurm",
}


#---------------------------------
# PBS STAGES
#---------------------------------
# The configuration options for individual stages. 

stages = {
    "makeIndexHtml": {
        "command": "python %html_index %name %samp %comp %results %output",
        "walltime": "10:00",
    },
    "fastQC": {
        "command": "fastqc -o %outdir -f fastq %pair1 %pair2",
        "walltime": "30:00"
    },
    "trimReads": {
        "command": "java -Xmx6g %tmp -jar %trimmomatic %paired -threads 1 " \
                   "-phred33 %log %parameters",
        "walltime": "10:00:00",
        "memInGB": 10,
    },
    "fastQCSummary": {
        "command": "python %script %fastqc_dir %fastqc_post_trim_dir " \
                   "%qc_summary %fastqc_summary %paired > %summary_txt",
        "walltime": "10:00",
    },
    "buildTranscriptomeIndex": {
        "command": "sh %buildIndex %seq %tmp_dir %index_dir %index %gene_ref " \
                   "%genome_ref",
        "walltime": "2:00:00",
        "queue": "smp"
    },
    "tophatAlign": {
        "command": "sh %tophat %pair1 %pair2 %out_dir %gene_ref %genome_ref " \
                   "%known %rgsm %rglb %rgid %rgpl %link",
        "walltime": "10:00:00",
        "memInGB": 32,
        "queue": "smp"
    },
    "sortBam": {
        "command": "samtools sort %bam %output",
        "walltime": "2:00:00"
    },
    "indexBam": {
        "command": "samtools index %bam",
        "walltime": "1:00:00"
    },
    "mergeTophat" : {
        "command": "sh %merge %fix_reads %samp_dir %output",
        "walltime": "1:00:00",
        "memInGB": 24
    },
    "reorderBam": {
        "command": "java -Xmx2g %tmp -jar %reorder_sam INPUT=%input " \
                   "OUTPUT=%output REFERENCE=%genome_ref",
        "walltime": "1:00:00"
    },
    "addRG": {
        "command": "java -Xmx2g %tmp -jar %add_rg INPUT=%bam OUTPUT=%output " \
                   "RGSM=%rgsm RGLB=%rglb RGID=%rgid RGPL=%rgpl RGPU=%rgpu",
        "walltime": "1:00:00"
    },
    "markDuplicates": {
        "command": "java -Xmx10g %tmp -jar %mark_dup INPUT=%input " \
                   "REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT " \
                   "AS=true METRICS_FILE=%log OUTPUT=%output",
        "walltime": "1:00:00",
        "memInGB": 12
    },
    "rnaSeQC": {
        "command": "java -Xmx10g %tmp -jar %rnaseqc %paired -n 1000 -s %samp " \
                   "-r %genome_ref -t %gene_ref %rrna_ref -o %outDir",
        "walltime": "1:00:00",
        "memInGB": 16
    },
    "cufflinksAssembly": {
        "command": "cufflinks -p 8 -o %outdir %in",
        "walltime": "2:00:00",
        "memInGB": 32,
        "queue": "smp"
    },
    "cuffmerge": {
        "command": "cuffmerge -g %gene_ref -s %genome_ref -p 8 -o %outdir %in",
        "walltime": "2:00:00",
        "memInGB": 32,
        "queue": "smp"
    },
    "cuffdiff": {
        "command": "cuffdiff %mask -o %outdir -p 8 -L %labels -u %merged_gtf " \
                   "%samples1 %samples2",
        "walltime": "8:00:00",
        "memInGB": 40,
        "queue": "smp"
    },
    "sortBamByName": {
        "command": "samtools sort -n %bam_file %sorted_file",
        "walltime": "1:00:00",
    },
    "alignmentStats": {
        "command": "sh %script %bam %unmapped %output %paired",
        "walltime": "1:00:00",
    },
    "qcSummary": {
        "command": "python %qc_script %fastqc_dir %fastqc_post_dir " \
                   "%alignment_stats_dir %rnaqc_dir %qc_summary %paired",
        "walltime": "10:00",
    },
    "countReads": {
        "command": "sh %htseq %bam %gene_ref %union %strict %stranded",
        "walltime": "5:00:00",
    },
    "combineAndAnnotate": {
        "command": "R --no-save --args %samples %comparisons %plain_text " \
                   "%rdata %annotate < %combine > %stdout 2> %stderr",
        "walltime": "30:00",
        "modules": ["R-intel/2.15.3"],
    },
    "voom": {
        "command": "R --no-save --args %rdata %outdir %voom < %script " \
                   "> %stdout 2> %stderr",
        "walltime": "30:00",
        "modules": ["R-intel/2.15.3"]
    },
    "edgeR": {
        "command": "R --no-save --args %rdata %outdir %edgeR < %script " \
                   "> %stdout 2> %stderr",
        "walltime": "30:00",
        "modules": ["R-intel/2.15.3"]
    },
}

