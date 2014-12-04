######################################################################
## Annotate genes using Ensembl BioMaRt
## Create CSV file of counts
## Set up factor of interest, covariates, and comparisons
######################################################################

args <- commandArgs(trailingOnly=TRUE)

sample.csv.filename <- args[1]
comparison.csv.filename <- args[2]
plain.text.counts <- args[3]
rdata.counts <- args[4]
biomart.dataset <- args[5]

print(args)

######################################################################
## Load data

samples <- read.table(sample.csv.filename, header=FALSE, 
                      stringsAsFactors=FALSE, sep=",",
                      colClasses="character")
comparisons <- read.table(comparison.csv.filename, header=FALSE,
                          stringsAsFactors=FALSE, sep=",",
                          colClasses="character")

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
## Load htseq count files

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
## Get annotations from Ensembl Biomart

if (biomart.dataset != "False") {
  library(biomaRt)
  mart <- useMart(biomart="ensembl", dataset=biomart.dataset)
  attribute.list <- listAttributes(mart)[,1]
  
  # Find which attribute name to get gene symbols from
  if ("hgnc_symbol" %in% attribute.list) {
    symbol <- "hgnc_symbol"
  } else if ("mgi_symbol" %in% attribute.list) {
    symbol <- "mgi_symbol"
  } else {
    # If symbol is unspecified, use the external_gene_id attribute or 
    # change it yourself by manually checking attribute list for your
    # organism's dataset
    symbol <- "external_gene_id"
  }
  sprintf("Using '%s' as symbol", symbol)
  
  # If you change the attributes, make sure there is only one entry per
  # ensembl_gene_id. The first two items should always be "ensembl_gene_id"
  # and symbol. If there is more than one entry per symbol, the first one 
  # is used (therefore don't include GO Terms in attributes).
  attributes <- c("ensembl_gene_id", symbol, "chromosome_name",
                  "start_position", "gene_biotype", "entrezgene", 
                  "description")
  print(attributes)
  annotations.raw <- getBM(attributes=attributes, mart=mart, verbose=T)
  
  # Remove duplicated entries
  annotations <- annotations.raw[!(duplicated(
                   annotations.raw$ensembl_gene_id)),]
  #duplicated <- annotations.raw[duplicated(annotations.raw$ensembl_gene_id),]
  
  gene.list <- rownames(counts)
  gene.list <- gene.list[gene.list %in% annotations$ensembl_gene_id]
  annotations <- annotations[annotations$ensembl_gene_id %in% gene.list,]
  rownames(annotations) <- annotations$ensembl_gene_id
  annotations <- annotations[,-1]
  
  # Only include counts which have an entry in ensembl
  counts <- counts[gene.list,]
  annotations <- annotations[gene.list,]
  
  # Check if annotations data frame is empty
  if (dim(annotations)[1] == 0) {
    # exit with error
    print("Error: annotation data frame is empty")
    quit("no", 1, FALSE)
  }
  if (any(rownames(counts) != rownames(annotations))) {
    print("Error: genes in counts and annotations don't match.")
    quit("no", 1, FALSE)
  }
  
} else {
  annotations <- c()
}


######################################################################
## Write to file

rm(i)
rm(args)

# Save to .RData file
save.image(rdata.counts)

# Save to txt file
# Add condition to sample names for columns in output file

counts.output <- cbind(Feature=rownames(counts), counts, annotations)
colnames(counts.output)[2:(dim(counts)[2]+1)] <- paste(sample.names,
                                                       condition, sep="_")
write.table(counts.output, plain.text.counts, col.names=TRUE, 
            row.names=FALSE, quote=FALSE, sep="\t")

