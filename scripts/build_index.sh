#!/bin/bash

# Run TopHat to create index file for reference transcriptome.

sample=$1
tmpDir=$2
indexDir=$3
index=$4
geneRef=$5
humanRef=$6

# Make mock fastq file to align
zcat $sample | head -1000 | gzip > ${tmpDir}/tmp.fastq.gz

# Make sure indexDir doesn't already exist. If it does, delete it.
if [ -d ${indexDir} ] && [ $(basename ${indexDir}) == "transcriptome_index" ]
    then
    rm -r ${indexDir}
fi

# Run TopHat and create transcriptome index
tophat \
    -p 8 \
    -G ${geneRef} \
    -o ${tmpDir} \
    --transcriptome-index=${index} \
    ${humanRef} \
    ${tmpDir}/tmp.fastq.gz


# Move log to indexDir and remove output tmpDir
mv ${tmpDir}/logs $indexDir

if [ $(basename $tmpDir) == "tmp" ]
    then 
    rm -r ${tmpDir}
fi


