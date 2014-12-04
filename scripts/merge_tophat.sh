#!/bin/bash

# Merge Tophat outputs accepted_hits.bam and unmapped.bam together for RNA-SeQC

fixTophatReads=$1
sampleDir=$2
mergedBam=$3

acceptedHits=${sampleDir}/accepted_hits.bam
unmapped=${sampleDir}/unmapped.bam
unmappedFixup=${sampleDir}/unmapped_fixup.bam

# Index accepted_hits.bam
samtools index $acceptedHits

# Fix unmapped reads (needs pysam library in Python)
python $fixTophatReads $sampleDir

# Merge accepted_hits.bam and unmapped_fixup.bam
samtools merge $mergedBam $acceptedHits $unmappedFixup

