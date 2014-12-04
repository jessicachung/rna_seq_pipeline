#!/bin/env python

# All possible columns:
# extract_fastqc = ["Conventional base calls", "Encoding", "Total Sequences", "Filtered Sequences", "Sequence length", "%GC"]
# extract_summary = ["Basic Statistics", "Per base sequence quality", "Per sequence quality scores", "Per base sequence content", "Per base GC content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels", "Overrepresented sequences", "Kmer Content"]

extract_fastqc = [
    "Total Sequences", 
    "%GC", 
    "Sequence length"]
extract_summary = [
    "Per base sequence quality", 
    "Per sequence GC content", 
    "Overrepresented sequences", 
    "Kmer Content"]

from sys import argv
from sys import exit
from glob import glob
import os.path

script, fastqc_dir, fastqc_post_trim_dir, fastqc_summary_file, \
        basic_statistics_file, paired_end = argv

paired_end = True if paired_end == "paired" else False

anchor_links = {"Basic Statistics": "M0",
                "Per base sequence quality": "M1",
                "Per sequence quality scores": "M2",
                "Per base sequence content": "M3",
                "Per base GC content": "M4",
                "Per sequence GC content": "M5",
                "Per base N content": "M6",
                "Sequence Length Distribution": "M7",
                "Sequence Duplication Levels": "M8",
                "Overrepresented sequences": "M9",
                "Kmer Content": "M10"}

CSS = """<html>
<head><title>FastQC Summary</title>
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

def parse_file(filename):
    file = open(filename)
    dict = {}
    for i in file.read().split("\n>>")[1:-1]:
        if i != "END_MODULE":
            lines = i.split("\n")
            module_name, status = lines[0].split("\t")
            dict[module_name] = lines
    file.close()
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

def parse_summary(module, extract):
    list = []
    for i in extract:
        try:
            list.append(module[i][0].split("\t")[1])
        except:
            list.append("-")
    return list

def print_list(list, extract, columns):
    print "\tPre-Trim" + "\t" * (len(extract_fastqc)-1) + "Post-Trim" + \
          "\t" * (len(extract_fastqc)-1)
    print "File\t" + "\t".join(extract * columns)
    for i in list:
        print "%s\t" % i[0],
        print "\t".join(i[1])
    print "\n"

def table_fastqc(output, list):
    output.write("<table>\n")
    output.write('<tr>\n<td></td>\n<td colspan="%d" align="center">' \
                 'FastQC Pre-Trim</td>\n<td colspan="%d" align="center">' \
                 'FastQC Post-Trim</td>\n</tr>' % (len(extract_summary),
                                                   len(extract_summary) ))
    output.write("<tr>\n<td>File</td>\n" + td(extract_summary, "left") * 2 + \
                 "</tr>\n" )
    for i in range(0,len(list)):
        output.write("<tr>\n")
        # file
        output.write(td([list[i][0]], "left"))
        # fastqc + hyperlinks
        output.write(td_pass_fail(list[i][1], "right", list[i][0]))
        output.write("</tr>\n")
    output.write("</table>")

def table_basic_statistics(output, list):
    output.write("<table>\n")
    output.write('<tr>\n<td></td>\n<td colspan="%d" align="center">' \
                 'FastQC Pre-Trim</td>\n<td colspan="%d" align="center">' \
                 'FastQC Post-Trim</td>\n' % (len(extract_fastqc),
                                              len(extract_fastqc)))
    output.write("<tr><td>File</td>\n" + td(extract_fastqc, "left") * 2 + \
                 "</tr>\n" )
    for i in range(0,len(list)):
        output.write("<tr>\n")
        # file
        output.write(td([list[i][0]], "left"))
        # fastqc
        output.write(td(list[i][1][0:(2*len(extract_fastqc))], "right"))
        output.write("</tr>\n")
    output.write("</table>")

def td(list, align):
    string = ""
    for i in list:
        string += "<td align=\"%s\">%s</td>\n" % (align, parse_number(i))
    return string

def td_pass_fail(list, align, sample_name):
    string = ""
    count = 0
    for i in list:
        if i == "pass":
            colour = "#C5D8A2"
        elif i == "warn":
            colour = "#FFFFE7"
        elif i == "fail":
            colour = "#FCD8D4"
        else:
            colour = "#FFFFFF"
        if count < len(extract_summary):
            link_dir = fastqc_dir
            suffix = "_fastqc"
        else:
            link_dir = fastqc_post_trim_dir
            if paired_end:
                suffix = ".trimmed-paired_fastqc"
            else:
                suffix = ".trimmed-single_fastqc"
        hyperlink = "%s/%s%s/fastqc_report.html" % (link_dir, sample_name, 
                                                    suffix)
        anchor = anchor_links[extract_summary[count % len(extract_summary)]]
        string += '<td align="%s" bgcolor="%s"><a href="%s#%s">%s</a>' \
                  '</td>\n' % (align, colour, hyperlink, anchor, i.upper())
        count += 1
    return string

def parse_number(number):
    try:
        int(number)
        return format(int(number), ",d")
    except:
        return number



    
def main():
    # Get fastQC files
    try:
        files = glob(os.path.join(fastqc_dir, "*/fastqc_data.txt"))
        post_trim_files = glob(os.path.join(fastqc_post_trim_dir, 
                               "*/fastqc_data.txt"))
    except:
        print "ERROR"
        exit()

    # Parse files
    samples = {}
    for filename in files:
        samples[filename] = parse_file(filename)
    post_trim_samples = {}
    for filename in post_trim_files:
        post_trim_samples[filename] = parse_file(filename)
        
        
    #---------------------------------------------
    # Parse module results: pass/warn/fail
    #---------------------------------------------
    module_results = {}
    for filename in samples:
        sample_name = filename.split("/")[-2][:-7]
        module_results[sample_name] = parse_summary(samples[filename], 
                                                    extract_summary)
    for filename in post_trim_samples:
        sample_name = filename.split("/")[-2][:-22]
        try:
            module_results[sample_name] = module_results[sample_name] + \
                    parse_summary(post_trim_samples[filename], extract_summary)
        except:
            pass

    # If no post-trim file, fill in empty cells with "-"
    for sample_name in module_results:
        if len(module_results[sample_name]) != 2 * len(extract_summary):
            module_results[sample_name] = module_results[sample_name] + \
                    ["-"] * len(extract_summary)

    # Print table to stdout
    module_sorted = module_results.items()
    module_sorted.sort()
    print "FastQC Summary"
    print_list(module_sorted, extract_summary, 2)
    
    # Output html table   
    try:
        output = open(fastqc_summary_file,'w')
        output.write(CSS)
        output.write("<body>\n<h1>FastQC Summary</h1>\n")
        table_fastqc(output, module_sorted)
        output.write("</body>\n</html>")
        output.close()
    except:
        print "ERROR. Could not create file %s." % fastqc_summary_file
    
    
    
    #----------------------------------------------------
    # Parse information from 'Basic Statistics' module
    #----------------------------------------------------
    basic_statistics_results = {}
    for filename in samples:
        sample_name = filename.split("/")[-2][:-7]
        basic_statistics_results[sample_name] = extract_info(
                samples[filename]["Basic Statistics"], extract_fastqc)
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
    
    # If no post-trim file, fill in empty cells with "-"
    for sample_name in basic_statistics_results:
        if len(basic_statistics_results[sample_name]) != \
                2 * len(extract_fastqc):
            basic_statistics_results[sample_name] = \
                    basic_statistics_results[sample_name] + \
                    ["-"] * len(extract_fastqc)

    
    # Print table to stdout
    basic_statistics_sorted = basic_statistics_results.items()
    basic_statistics_sorted.sort()
    print "FastQC Basic Statistics Summary"
    print_list(basic_statistics_sorted, extract_fastqc, 2)
    
    # Output html table
    try:
        output = open(basic_statistics_file,'w')
        output.write(CSS)
        output.write("<body>\n<h1>FastQC Basic Statistics Summary</h1>\n")
        table_basic_statistics(output, basic_statistics_sorted)
        output.write("</body>\n</html>")
        output.close()
    except:
        print "ERROR. Could not create file %s." % basic_statistics_file
    

if __name__ == "__main__":
    main()


