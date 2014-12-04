#!/bin/bash

# Count gene features using HTSeq-count

bamFile=$1
sortedSamFile=${bamFile:0:${#bamFile}-4}.sam
geneRef=$2
htseqCount=$3
htseqStrict=$4
stranded=$5


samtools view -h ${bamFile} > ${sortedSamFile}

# union HTSeq
htseq-count \
    --stranded=${stranded} \
    ${sortedSamFile} \
    ${geneRef} \
    > ${htseqCount}

# intesection-strict HTSeq
htseq-count \
    -m intersection-strict \
    --stranded=${stranded} \
    ${sortedSamFile} \
    ${geneRef} \
    > ${htseqStrict}

rm ${sortedSamFile}

