#!/bin/env python

# All possible columns:
'''
extract_fastqc = [
    "Conventional base calls", 
    "Encoding", 
    "Total Sequences", 
    "Filtered Sequences", 
    "Sequence length", 
    "%GC"]
extract_alignment = [
    "Number of left reads", 
    "Number of right reads", 
    "Paired reads with only one aligned pair", 
    "Paired reads with unique alignments", 
    "Paired reads with multiple alignments", 
    "Paired reads with only one unaligned pair", 
    "Paired reads which could not be aligned"]
extract_rnaseqc = [
    "Sample", 
    "Note", 
    "End 2 Mapping Rate", 
    "Chimeric Pairs", 
    "Intragenic Rate", 
    "Num. Gaps", 
    "Mapping Rate", 
    "Exonic Rate", 
    "5' Norm", 
    "Genes Detected", 
    "Unique Rate of Mapped", 
    "Read Length", 
    "Mean Per Base Cov.", 
    "End 1 Mismatch Rate", 
    "Fragment Length StdDev", 
    "Estimated Library Size", 
    "Mapped", 
    "Intergenic Rate", 
    "rRNA", 
    "Total Purity Filtered Reads Sequenced", 
    "Failed Vendor QC Check", 
    "Mean CV", 
    "Transcripts Detected", 
    "Mapped Pairs", 
    "Cumul. Gap Length", 
    "Gap %", 
    "Unpaired Reads", 
    "Intronic Rate", 
    "Mapped Unique Rate of Total", 
    "Expression Profiling Efficiency", 
    "Mapped Unique", 
    "End 2 Mismatch Rate", 
    "End 2 Antisense", 
    "Alternative Aligments", 
    "End 2 Sense", 
    "Fragment Length Mean", 
    "End 1 Antisense", 
    "Base Mismatch Rate", 
    "End 1 Sense", 
    "End 1 % Sense", 
    "rRNA rate", 
    "End 1 Mapping Rate", 
    "No. Covered 5'", 
    "Duplication Rate of Mapped", 
    "End 2 % Sense"]
'''

from sys import argv
from sys import exit
from glob import glob
import os.path
import re

script, fastqc_dir, fastqc_post_trim_dir, alignment_stats_dir, rnaseqc_dir, \
        qc_summary_file, paired_end = argv


extract_fastqc = ["Total Sequences"]

if paired_end == "paired":
    extract_alignment = [
        "Number of left reads", 
        "Number of right reads", 
        "Paired reads with unique alignments", 
        "Paired reads with multiple alignments", 
        "Paired reads which could not be aligned"]
else:
    extract_alignment = [
        "Number of reads", 
        "Reads with unique alignments", 
        "Reads with multiple alignments", 
        "Reads which could not be aligned"]

extract_rnaseqc = [
    "Mapped", 
    "Mapping Rate",
    "Mapped Pairs",
    "Unpaired Reads",
    "Intragenic Rate",
    "Exonic Rate",
    "Mapped Unique",
    "Unique Rate of Mapped",
    "Mapped Unique Rate of Total",
    "Duplication Rate of Mapped"]

CSS = """<html>
<head><title>RNA-Seq QC Summary</title>
<style type="text/css">
table {
        border-width: 1px;
        border-spacing: 2px;
        border-style: solid;
        border-color: gray;
        border-collapse: collapse;
}
table td {
        border-width: 2px;
        padding: 4px;
        border-style: solid;
        border-color: gray;
}
</style>
</head>
"""

def parse_fastqc(filename):
    file = open(filename)
    dict = {}
    for i in file.read().split("\n>>")[1:-1]:
        if i != "END_MODULE":
            lines = i.split("\n")
            module_name, status = lines[0].split("\t")
            dict[module_name] = lines
    file.close()
    return dict

def parse_alignment_stats(filename):
    file = open(filename)
    dict = {}
    for line in file.read().strip().split("\n"):
        metric = line.split("\t")
        dict[metric[1]] = metric[0]
    file.close()
    return dict

def parse_rna_file(filename):
    file = open(filename)
    dict = {}
    lines = file.read().split("\n")
    metrics = lines[0].split("\t")
    values = lines[1].split("\t")
    for i in range(0,len(metrics)):
        dict[metrics[i]] = values[i]
    return dict

def extract_info(module, extract):
    dict = {}
    list = []
    for i in module:
        dict[i.split("\t")[0]] = i.split("\t")[1]
    for i in extract:
        try:
            list.append(dict[i])
        except:
            list.append("-")
    return list

def extract_alignment_stats(dict, extract):
    list = []
    for i in extract:
        try:
            list.append(dict[i])
        except:
            list.append("-")
    return list

def table(output, basic_statistics, list):
    output.write("<table>\n")
    output.write('<tr>\n<td colspan=2></td>\n<td colspan="%d" ' \
                 'align="center">FastQC Pre-Trim</td>\n<td colspan="%d" ' \
                 'align="center">FastQC Post-Trim</td>\n<td colspan="%d" ' \
                 'align="center">TopHat Alignment Stats</td>\n<td ' \
                 'colspan="%d" align="center">RNA-SeQC</td>\n</tr>' % \
                 (len(extract_fastqc), len(extract_fastqc),
                  len(extract_alignment), len(extract_rnaseqc)))
    output.write("<tr><td>Sample</td>\n<td>File</td>\n" + \
                 td(extract_fastqc, "left") * 2 + \
                 td(extract_alignment, "left") + \
                 td(extract_rnaseqc, "left") + "</tr>\n")
    for i in range(0,len(list)):
        row_span = len(list[i][1][0])
        output.write("<tr>\n")
        # Sample
        output.write(td_rowspan([list[i][0]], row_span, "left"))
        count = 0
        for file in list[i][1][0]:
            if count == 0:
                # Files
                output.write(td([file], "left"))
                # FastQC stats
                output.write(td(basic_statistics[file], "right"))
                # Alignment and RNA-SeQC stats
                output.write(td_rowspan(list[i][1][1], row_span, "right"))
            else:
                output.write("<tr>\n")
                # Files
                output.write(td([file], "left"))
                # FastQC stats
                output.write(td(basic_statistics[file], "right"))
            count += 1
            output.write("</tr>\n")
    output.write("</table>")

def td(list, align):
    string = ""
    for i in list:
        string += '<td align="%s">%s</td>\n' % (align, parse_number(i))
    return string

def td_colour(list, align):
    string = ""
    for i in list:
        if i == "pass":
            colour = "#C5D8A2"
        elif i == "warn":
            colour = "#FFFFE7"
        elif i == "fail":
            colour = "#FCD8D4"
        else:
            colour = "#FFFFFF"
        string += '<td align="%s" bgcolor="%s">%s</td>\n' % \
                (align, colour, i.upper())
    return string

def td_rowspan(list, row_span, align):
    string = ""
    if row_span == None:
        return string
    for i in range(0,len(list)):
        try:
            string += '<td rowspan=%d align="%s">%s</td>\n' % \
                    (row_span, align, parse_number(list[i]))
        except:
            string += '<td rowspan=%d align="%s">%s</td>\n' % \
                    (row_span, align, "-")
    return string

def parse_number(number):
    try:
        int(number)
        return format(int(number), ",d")
    except:
        return number



    
def main():

    try:
        pre_trim_files = glob(os.path.join(fastqc_dir, "*/fastqc_data.txt"))
        post_trim_files = glob(os.path.join(fastqc_post_trim_dir, 
                                            "*/fastqc_data.txt"))
        alignment_stats_files = glob(os.path.join(alignment_stats_dir,
                                                  "*.alignmentStats.txt"))
        rnaqc_files = glob(os.path.join(rnaseqc_dir, "*/metrics.tsv"))
    except:
        print "ERROR"
        exit()

    #-------------------------------------------------------
    # FastQC stats for pre-trimming and post-trimming
    #-------------------------------------------------------
    # parse files
    pre_trim_samples = {}
    for filename in pre_trim_files:
        pre_trim_samples[filename] = parse_fastqc(filename)
    post_trim_samples = {}
    for filename in post_trim_files:
        post_trim_samples[filename] = parse_fastqc(filename)
    
    # extract info
    basic_statistics_results = {}
    for filename in pre_trim_samples:
        sample_name = filename.split("/")[-2][:-7]
        basic_statistics_results[sample_name] = extract_info(
                pre_trim_samples[filename]["Basic Statistics"], extract_fastqc)
    information_post_trim = {}
    for filename in post_trim_samples:    
        sample_name = filename.split("/")[-2][:-22]
        try:
            basic_statistics_results[sample_name] = \
                    basic_statistics_results[sample_name] + \
                    extract_info(post_trim_samples[filename]\
                                 ["Basic Statistics"], extract_fastqc)
        except:
            pass
    
    # if no post-trim file, fill in empty cells with "-"
    for sample_name in basic_statistics_results:
        if len(basic_statistics_results[sample_name]) != \
                2 * len(extract_fastqc):
            basic_statistics_results[sample_name] = \
                    basic_statistics_results[sample_name] + \
                    ["-"] * len(extract_fastqc)


    #-------------------------------------------------------
    # Alignment stats from TopHat
    #-------------------------------------------------------
    # parse files
    alignment_samples = {}
    for filename in alignment_stats_files:
        sample_name = filename.split("/")[-1][:-19]
        alignment_samples[sample_name] = parse_alignment_stats(filename)
    
    # extract info
    alignment_stats_results = {}
    for sample_name in alignment_samples:
        alignment_stats_results[sample_name] = \
                extract_alignment_stats(alignment_samples[sample_name], 
                                        extract_alignment)


    #-------------------------------------------------------
    # RNA-Seq stats from RNA-SeQC
    #-------------------------------------------------------
    # parse RNA-SeQC files
    rna_samples = {}
    for filename in rnaqc_files:
        sample_name = filename.split("/")[-2]
        rna_samples[sample_name] = parse_rna_file(filename)

    # extract info
    for sample_name in rna_samples:
        alignment_stats_results[sample_name] = \
                alignment_stats_results[sample_name] + \
                extract_alignment_stats(rna_samples[sample_name], 
                                        extract_rnaseqc)
    

    #-------------------------------------------------------
    # Output to HTML table
    #-------------------------------------------------------    
    # join dictionaries
    statistics = {}
    for sample_name in alignment_stats_results:
        individual_files = []
        for file in basic_statistics_results:
            match = re.search('(SM_[A-Za-z0-9-.]+_RP_[A-Za-z0-9-.]+)_.*', 
                              os.path.basename(file)).group(1)
            if sample_name == match:
                individual_files.append(file)
        individual_files.sort()
        statistics[sample_name] = [individual_files,
                                   alignment_stats_results[sample_name]]
  
    statistics_sorted = statistics.items()
    statistics_sorted.sort()
    
    try:
        output = open(qc_summary_file,'w')
        output.write(CSS)
        output.write("<body>\n<h1>QC Metrics</h1>\n")
        table(output, basic_statistics_results, statistics_sorted)
        output.write("</body>\n</html>")
        output.close()
    except:
        print "ERROR. Could not create file %s." % qc_summary_file
    
    


if __name__ == "__main__":
    main()

