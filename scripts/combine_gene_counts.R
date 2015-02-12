######################################################################
## Combine counts from HTSeq into one matrix
## Annotate genes using Ensembl BioMaRt
## Create CSV file of counts
## Set up factor of interest, covariates, and comparisons
######################################################################

args <- commandArgs(trailingOnly=TRUE)

sample.csv.filename <- args[1]
comparison.csv.filename <- args[2]
plain.text.counts <- args[3]
rdata.counts <- args[4]
annotations.data <- args[5]

print(args)

######################################################################
## Load data

samples <- read.table(sample.csv.filename, header=FALSE, 
                      stringsAsFactors=FALSE, sep=",",
                      colClasses="character")
comparisons <- read.table(comparison.csv.filename, header=FALSE,
                          stringsAsFactors=FALSE, sep=",",
                          colClasses="character")
load(annotations.data)

# Check if annotations available
if (length(annotations) == 0) {
  annotate <- FALSE
} else {
  annotate <- TRUE
}

######################################################################
## Set up conditions and covariates

sample.names <- samples[,1]
htseq.filenames <- samples[,2]

# Condition of interest
condition <- factor(samples[,3])

n.covariates <- ncol(samples) - 3
n.comparisons <- nrow(comparisons)

# Covariates
if (n.covariates > 0) {
  covariates <- samples[,4:ncol(samples),drop=F]
} else{
  covariates <- FALSE
}

######################################################################
## Load HTSeq count files

counts <- c()
for (i in samples[,2]){
  counts <- cbind(counts, read.table(i, header=FALSE)[,2])
}
rownames(counts) <- read.table(samples[1,2], header=FALSE)[,1]
colnames(counts) <- samples[,1]

# Remove last rows from htseq-count counts
counts <- counts[1:(dim(counts)[1]-5),]

head(counts)
tail(counts)

# Remove genes with no counts (comment out if you want all genes)
counts <- counts[rowSums(counts) != 0,]


######################################################################
## Annotate counts matrix with Ensembl data

if (annotate) {
  gene.list <- rownames(counts)
  gene.list <- gene.list[gene.list %in% rownames(annotations)]
  
  # Only include counts which have an entry in ensembl
  counts <- counts[gene.list,]
  annotations <- annotations[gene.list,]
}


######################################################################
## Write to file

# Save to .RData file
save.image(rdata.counts)

# Save to txt file
# Add condition to sample names for columns in output file
counts.output <- cbind(Feature=rownames(counts), counts, annotations)
colnames(counts.output)[2:(dim(counts)[2]+1)] <- paste(sample.names,
                                                       condition, sep="_")
write.table(counts.output, plain.text.counts, col.names=TRUE, 
            row.names=FALSE, quote=FALSE, sep="\t")


proc.time()
sessionInfo()
