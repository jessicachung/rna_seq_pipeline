#!/usr/bin/env python
# make index.html in results directory

import sys
from glob import glob
import os.path
from copy import deepcopy
import re

script, analysis_name, sample_csv_file, comparison_csv_file, output_dir, \
        output_filename = sys.argv


# Output directories:
fastqc_dir = os.path.join(output_dir, "fastqc")
fastqc_post_trim_dir = os.path.join(output_dir, "fastqc_post_trim")
trimmed_dir = os.path.join(output_dir, "trimmed_reads")
tophat_raw_dir = os.path.join(output_dir, "tophat_raw")
tophat_dir = os.path.join(output_dir, "tophat")
cufflinks_dir = os.path.join(output_dir, "cufflinks")
cuffmerge_dir = os.path.join(output_dir, "cuffmerge")
cuffdiff_dir = os.path.join(output_dir, "cuffdiff")
htseq_dir = os.path.join(output_dir, "htseq_count")
counts_dir = os.path.join(output_dir, "read_counts")
merged_dir = os.path.join(output_dir, "tophat_merged")
rnaseqc_dir = os.path.join(output_dir, "rnaseqc")
alignment_stats_dir = os.path.join(output_dir, "alignment_stats")
qc_summary_dir = os.path.join(output_dir, "qc_summary")

main_voom_dir = os.path.join(output_dir, "voom_analysis")
main_edger_dir = os.path.join(output_dir, "edgeR_analysis")
voom_dir = os.path.join(main_voom_dir, analysis_name + "_voom")
edger_dir = os.path.join(main_edger_dir, analysis_name + "_edgeR")
cuffmerge_sub_dir = os.path.join(cuffmerge_dir, analysis_name + "_cuffmerge")


##########################################################################


def html_table(heading, list_table):
    html_output = []
    html_output.append("<hr width=80%>")
    html_output.append("<h2>%s</h2>" % heading)
    html_output.append("<table width=90%>")
    for list in list_table:
        html_output.append("<tr>")
        html_output.append("<td width=20%%><a href=\"%s\">%s</a></td>" % \
                           (list[0],list[1]))
        html_output.append("<td>%s</td>" % list[2])
        html_output.append("</tr>")
    html_output.append("</table>")
    return "\n".join(html_output)

def no_output(heading, message):
    html_output = []
    html_output.append("<hr width=80%>")
    html_output.append("<h2>%s</h2>" % heading)
    html_output.append("<p>%s<p>" % message)
    return "\n".join(html_output)


##########################################################################
## Summary


CSS = """<html>
<head><title>RNA-Seq Results Index</title>
<style type="text/css">
table {
    border-width: 1px;
    border-spacing: 2px;
    border-style: solid;
    border-color: gray;
    border-collapse: collapse;
    font: 12pt Georgia;
}
table td {
    border-width: 2px;
    padding: 4px;
    border-style: solid;
    border-color: gray;
}
p.summary {
	background-color: #FAF2D2;
	padding: 10px 10px 10px 10px;
	border: 1px #FFD54D solid;
	width: 80%;
	text-align: left !important;
}
body {
    font: 10pt Verdana;
	margin: 10px 100px;
}
</style>
</head>
<body>
<center>"""

summary_html = """
<h2>Index</h2>
<p class=summary>
This HTML file was generated with the sample CSV file:
<br>
<a href=\"%s\">%s</a>
<br>
and the comparison CSV file:
<br>
<a href=\"%s\">%s</a>
<br>
and the results directory:
<br>
<a href=\"%s\">%s</a>
<br><br>
To update this file, rerun the 'makeIndexHtml' stage:
<br>
eg.
<br>
<code>
python RNA-seq_pipeline.py \\
<br>--opts=pipeline_config,pipeline_stages_config \\
<br>--style=run \\
<br>--end=makeIndexHtml \\
<br>--force=makeIndexHtml
</code>
</p>
""" % (sample_csv_file, sample_csv_file, comparison_csv_file, 
       comparison_csv_file, output_dir, output_dir)


##########################################################################
## QC

# fastQC pre and post trim summary, qc summary file

qc_files = glob(qc_summary_dir + "/*")
qc_files = map(lambda x: os.path.basename(x), qc_files)
qc_files_dict = {
    0: ["FastQC_basic_statistics_summary.html", 
        "FastQC Basic Statistics Summary", 
        "Statistics for each FASTQ file (pre and post-trimming) including read counts, %GC content, and sequence length."],
    1: ["FastQC_summary.html", 
        "FastQC Summary", 
        "Pass/Fail summary for all FASTQ files (pre and post-trimming) for metrics such as sequence quality, GC content, and overrepresented sequences"],
    2: ["qc_summary.html", 
        "RNA-Seq QC Summary", 
        "QC statistics for each sample. Includes TopHat alignment statistics and RNASeQC statistics."],
    3: ["rnaseqc/", 
        "RNASeQC Directory", 
        "Detailed RNASeQC statistics for each sample."]
    }

qc_table = []
relative_path = os.path.basename(qc_summary_dir) + "/"
for i in range(0,len(qc_files_dict)):
    if qc_files_dict[i][0] in qc_files:
        qc_files_dict[i][0] = relative_path + qc_files_dict[i][0]
        qc_table.append(qc_files_dict[i])


rnaseqc_success_files = glob(rnaseqc_dir + "/*.rnaSeQC.Success")

# Individual RNASeQC results
if len(rnaseqc_success_files) > 0:
    for i in range(0,len(qc_files_dict)):
        if qc_files_dict[i][0] == "rnaseqc/":
            qc_table.append(qc_files_dict[i])

if len(qc_table) > 0:
    qc_html = html_table("QC Results", qc_table)
else:
    qc_html = no_output("QC Results", "No QC output available.")


##########################################################################
## Read counts

read_count_files = glob(counts_dir + "/*")
read_count_files = map(lambda x: os.path.basename(x), read_count_files)
samples_csv_name = analysis_name + "_counts.txt"

if samples_csv_name in read_count_files:
    read_counts_table = [[os.path.join(counts_dir,samples_csv_name), 
                          "Read counts", 
                          "Read counts in plain text format. Includes annotations if provided. Can be uploaded for analysis in <a href=\"http://www.vicbioinformatics.com/degust/\">Degust</a>."]]
else:
    read_counts_table = []

if len(read_counts_table) > 0:
    read_counts_html = html_table("Read Counts", read_counts_table)
else:
    read_counts_html = no_output("Read Counts", "No counts available.")


##########################################################################
## Voom results

voom_files = glob(voom_dir + "/*")
voom_files = map(lambda x: os.path.basename(x), voom_files)
voom_files_dict = {
    0: ["top_genes_.*.html", 
        "Top 100 DE Genes (HTML)", 
        "Top 100 differentially expressed genes from Voom in HTML table format. Includes raw read counts for each gene from each sample (and TMM-normalised counts per million reads in parentheses). Includes gene annotations if provided."], 
    1: ["top_genes_.*.txt", 
        "Top 100 DE Genes", 
        "Top 100 differentially expressed genes from Voom in plain text format. Includes gene annotations if provided."],
    2: ["all_genes_.*.txt", 
        "All Genes", 
        "All genes from Voom sorted by descending p-value in plain text format."],
    3: ["plots_.*.pdf", 
        "DE Plots", 
        "Page 1: MA plot.\nPage 2: P-value distribution\nPage 3 & 4: Heatmap\nPage 5+: Histograms"],
    4: ["plots.pdf", 
        "QC Plots", 
        "Page 1: Raw gene count boxplots.\nPage 2: TMM-normalised gene count boxplots.\nPage 3: RLE plot of raw counts.\nPage 4: RLE plot of TMM-normalised counts.\nPage5: Voom mean-variance trend plot."],
    5: ["MDS.pdf", 
        "MDS Plots", 
        "MDS plots using multiple dimensions."],
    6: ["PCA.pdf", 
        "PCA Plots", 
        "PCA pair plot with 4 dimensions."],
    7: ["heatmap.pdf", 
        "Clustered Heatmap", 
        "Heatmap of the top 1000 genes with most variance."],
    8: ["voom.stdout", 
        "Voom Stdout file", 
        "Standard output from R."]
    }

voom_table = []
relative_path = "/".join(voom_dir.split("/")[-2:]) + "/"
if "voom.Success" in voom_files:
    for i in range(0,len(voom_files_dict)):
        for j in range(0,len(voom_files)):
            if re.search(voom_files_dict[i][0], voom_files[j]):
                new_row = deepcopy(voom_files_dict[i])
                new_row[0] = relative_path + voom_files[j]
                new_row[1] = voom_files[j]
                voom_table.append(new_row)
    voom_html = html_table("Voom Results", voom_table)
else:
    voom_html = no_output("Voom Results", "No output from Voom available.")


##########################################################################
## edgeR results

edger_files = glob(edger_dir + "/*")
edger_files = map(lambda x: os.path.basename(x), edger_files)
edger_files_dict = {
    0: ["top_genes_.*.html", 
        "Top 100 DE Genes (HTML)", 
        "Top 100 differentially expressed genes from edgeR in HTML table format. Includes raw read counts for each gene from each sample (and TMM-normalised counts per million reads in parentheses). Includes gene annotations if provided."], 
    1: ["top_genes_.*.txt", 
        "Top 100 DE Genes", 
        "Top 100 differentially expressed genes from edgeR in plain text format. Includes gene annotations if provided."],
    2: ["all_genes_.*.txt", 
        "All Genes", 
        "All genes from edgeR sorted by descending p-value in plain text format."],
    3: ["plots_.*.pdf", 
        "DE Plots", 
        "Page 1: MA plot.\nPage 2: P-value distribution\nPage 3 & 4: Heatmap\nPage 5+: Histograms"],
    4: ["plots.pdf", 
        "QC Plots", 
        "Page 1: Raw gene count boxplots.\nPage 2: TMM-normalised gene count boxplots.\nPage 3: RLE plot of raw counts.\nPage 4: RLE plot of TMM-normalised counts.\nPage5: edgeR BCV plot."],
    5: ["MDS.pdf", 
        "MDS Plots", 
        "MDS plots using multiple dimensions."],
    6: ["PCA.pdf", 
        "PCA Plots", 
        "PCA pair plot with 4 dimensions."],
    7: ["heatmap.pdf", 
        "Clustered Heatmap", 
        "Heatmap of the top 1000 genes with most variance."],
    8: ["edgeR.stdout", 
        "edgeR Stdout file", 
        "Standard output from R."]
    }

edger_table = []
relative_path = "/".join(edger_dir.split("/")[-2:]) + "/"
if "edgeR.Success" in edger_files:
    for i in range(0,len(edger_files_dict)):
        for j in range(0,len(edger_files)):
            if re.search(edger_files_dict[i][0], edger_files[j]):
                new_row = deepcopy(edger_files_dict[i])
                new_row[0] = relative_path + edger_files[j]
                new_row[1] = edger_files[j]
                edger_table.append(new_row)
    edger_html = html_table("EdgeR Results", edger_table)
else:
    edger_html = no_output("EdgeR Results", "No output from EdgeR available.")


##########################################################################
## Cuffdiff results

cuffdiff_success_files = glob(cuffdiff_dir + "/" + analysis_name + 
                              "*.cuffdiff.Success")
cuffdiff_sub_dirs = glob(cuffdiff_dir + "/" + analysis_name + "_*/")
cuffdiff_files_dict = {
    0: ["gene_exp.diff", 
        "Cuffdiff Expression Results", 
        "Gene expression results from Cuffdiff"]
    }

if len(cuffdiff_success_files) > 0:
    cuffdiff_table = []
    for dir in cuffdiff_sub_dirs:
        cuffdiff_files = glob(dir + "/*")
        cuffdiff_files = map(lambda x: os.path.basename(x), cuffdiff_files)
        for i in range(0,len(cuffdiff_files_dict)):
            if cuffdiff_files_dict[i][0] in cuffdiff_files:
                new = deepcopy(cuffdiff_files_dict[i])
                new[0] = "%s/%s" % (dir, cuffdiff_files_dict[i][0])
                new[2] += " (%s)." % os.path.basename(dir[:-1])
                cuffdiff_table.append(new)
    cuffdiff_html = html_table("Cuffdiff Results", cuffdiff_table)
else:
    cuffdiff_html = no_output("Cuffdiff Results", 
                              "No output from Cuffdiff available.")



##########################################################################
## Other directories

misc_dir = glob(output_dir + "/*")
misc_dir = map(lambda x: os.path.basename(x), misc_dir)

misc_dir_dict = {
    0: ["trimmed_reads", 
        "Trimmed Reads Directory", 
        "Contains FASTQ files with adapters trimmed by Trimmomatic."],
    1: ["fastqc", 
        "Pre-trimmed FastQC Directory", 
        "Contains FastQC results for FASTQ files before trimming."],
    2: ["fastqc_post_trim", 
        "Post-trimmed FastQC Directory", 
        "Contains FastQC results for FASTQ files after trimming."],
    3: ["transcriptome_index", 
        "Transcriptome Index Directory", 
        "Contains transcriptome index files built from reference files for TopHat to use for alignment."],
    4: ["tophat_raw", 
        "Raw TopHat Output Directory", 
        "Contains TopHat output."],
    5: ["tophat", 
        "Tophat Directory", 
        "Contains sorted and indexed TopHat alignment files."],
    6: ["tophat_merged", 
        "Merged TopHat Directroy", 
        "Contains TopHat output with accepted hits and unmapped reads merged into one alignment file. Also contains outputs from processes such as reordering, adding read groups, and duplicate flaging which are needed for RNASeQC."],
    7: ["alignment_stats", 
        "Alignment Stats Directory", 
        "Contains alignment statistics of TopHat alignment files."],
    8: ["rnaseqc", 
        "RNASeQC Directory", 
        "Contains RNASeQC output."],
    9: ["cufflinks", 
        "Cufflinks Directory", 
        "Contains Cufflinks assemblies."],
    10: ["cuffmerge",
         "Cuffmerge Directory",
         "Contains merged assemblies by Cuffmerge."],
    11: ["cuffdiff",
         "Cuffdiff Directory",
         "Contains differential gene expression analysis by Cuffdiff."],
    12: ["htseq_count",
         "HTSeq Directory",
         "Contains read counts for each gene from HTSeq-count."],
    13: ["read_counts",
         "Read Counts Directory",
         "Contains read counts filtered and merged into one file. Includes gene annotations if available."],
    14: ["voom_analysis",
         "Voom Directory",
         "Contains output from Voom analysis."],
    15: ["edgeR_analysis",
         "edgeR Directory",
         "Contains output from edgeR analysis."],
    }
    

misc_table = []
for i in range(0,len(misc_dir_dict)):
    if misc_dir_dict[i][0] in misc_dir:
        misc_table.append(misc_dir_dict[i])

misc_html = html_table("Result Directories", misc_table)



##########################################################################
## Output

html = [CSS, summary_html, qc_html, read_counts_html, voom_html, edger_html, 
        cuffdiff_html, misc_html, "</center></body>\n</html>"]


output_file = open(output_filename, 'w')
for i in html:
    output_file.write(i)
    output_file.write("\n<br><br>\n")

output_file.close()




