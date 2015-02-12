######################################################################
## Get annotations from Ensembl BioMaRt
######################################################################

args <- commandArgs(trailingOnly=TRUE)

biomart.dataset <- args[1]
biomart.dataset

output.dir <- args[2]
output.dir


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
  annotations.raw <- getBM(attributes=attributes, mart=mart, verbose=FALSE)
  
  # Remove duplicated entries
  annotations <- annotations.raw[!(duplicated(
                   annotations.raw$ensembl_gene_id)),]
  #duplicated <- annotations.raw[duplicated(annotations.raw$ensembl_gene_id),]
  
  rownames(annotations) <- annotations$ensembl_gene_id
  annotations <- annotations[,-1]
    
  # Check if annotations data frame is empty
  if (dim(annotations)[1] == 0) {
    # exit with error
    print("Error: annotation data frame is empty")
    quit("no", 1, FALSE)
  }
  
} else {
  annotations <- c()
}


######################################################################
## Write to file

# Save to .RData file
rdata.output <- paste(output.dir, "/annotations.RData", sep="")
save.image(rdata.output)


proc.time()
sessionInfo()
