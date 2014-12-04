#!/bin/bash

# Resume unfinished TopHat then create symbolic link for accepted_hits.bam

pair1=$1
pair2=$2
outDir=$3
geneRef=$4
humanRef=$5
known=$6
rgsm=$7
rglb=$8
rgid=$9
rgpl=${10}
symLink=${11}
outBam=${outDir}/accepted_hits.bam

# Run TopHat
tophat --resume ${outDir}


# Create relative symlink to accepted_hits.bam
tophatDir=`dirname $symLink`
relativeLink=`echo $outDir | awk -F "/" '{print "../" $(NF-1) "/" $NF "/accepted_hits.bam"}'`
newLinkName=`basename ${symLink}`

cd ${tophatDir}
ln -sf ${relativeLink} ${newLinkName}


# # Create absolute link to accepted_hits.bam
# outBam=${outDir}/accepted_hits.bam
# ln -sf ${outBam} ${symLink}
