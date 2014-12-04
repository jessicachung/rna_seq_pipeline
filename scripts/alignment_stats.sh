#!/bin/bash

# Get alignment stats including number of reads, number of reads which align
# uniquely, and number of reads with multiple alignments.

bamFile=$1
unmappedBam=$2
readStats=`dirname $2`/prep_reads.info
output=$3
paired=$4

if [[ ${paired} == "paired" ]]
then
    leftReads=`grep left_reads_out ${readStats} | cut -d "=" -f 2`
    rightReads=`grep right_reads_out ${readStats}    | cut -d "=" -f 2`
    echo -e "${leftReads}\tNumber of left reads" > ${output}
    echo -e "${rightReads}\tNumber of right reads" >> ${output}

    samtools view ${bamFile} \
        | cut -f 1 \
        | uniq -c \
        | awk 'BEGIN {x=0;y=0;z=0} 
               {if ($1 == 1) {x++} else if ($1 == 2) {y++} else {z++}}
               END {print x "\tPaired reads with only one aligned pair\n" \
                    y "\tPaired reads with unique alignments\n" \
                    z "\tPaired reads with multiple alignments"}' \
        >> ${output}

    samtools view ${unmappedBam} \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | awk 'BEGIN {x=0;y=0} 
               { if ($1 == 1) {x++} else {y++} }
               END {print x "\tPaired reads with only one unaligned pair\n" \
                    y "\tPaired reads which could not be aligned"}' \
        >> ${output}
else
    reads=`grep reads_out $readStats | cut -d "=" -f 2`
    echo -e "${reads}\tNumber of reads" > ${output}
    
    samtools view ${bamFile} \
        | cut -f 1 \
        | uniq -c \
        | awk 'BEGIN {x=0;y=0} 
               { if ($1 == 1) {x++} else {y++} }
               END {print x "\tReads with unique alignments\n" \
                    y "\tReads with multiple alignments"}' \
        >> ${output}
    
    samtools view ${unmappedBam} \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | awk 'BEGIN {x=0;y=0} 
               { if ($1 == 1) {x++} else {y++} }
               END {print x "\tReads which could not be aligned"}' \
        >> ${output}
fi
