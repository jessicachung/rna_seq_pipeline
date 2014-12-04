######################################################################
## Differential expression analysis script for pipeline
## Specify either voom or edgeR
## Fit data using a linear model to find differentially expressed genes.
## Design matrix is [ ~0+condition ] if no covariates 
##   or [ ~0+condition+c1+...+cn ] for n covariates.
## Uses a comparison matrix to make comparisons
######################################################################

library(limma)
library(edgeR)
library(statmod)
library(reshape)
library(ggplot2)
library(scales)
library(gridExtra)
library(xtable)

args <- commandArgs(trailingOnly=TRUE)

count.rdata <- args[1]
count.rdata

output.dir <- args[2]
output.dir

mode <- args[3]
mode

voom <- ifelse(mode == "voom", TRUE, FALSE)

# Number of genes to output to top_genes_c1_vs_c2.txt with annotations. Set
# to Inf to output all genes.
n.top <- 100

setwd(output.dir)
options(digits=3)

######################################################################
## Functions

VoomPlotMA <- function(fit, comparison) {
  values <- topTable(fit, coef=comparison, adjust="BH", sort="p", n=Inf)

  # Center plots around zero  
  max.M <- max(abs(min(values$logFC)), abs(max(values$logFC)))
  lim <- c(max.M * -1, max.M)
  
  # Highlight p < 0.05 genes in red
  col <- ifelse(values$adj.P.Val < 0.05, "red", "black")
  
  plot(values$AveExpr, values$logFC, pch=".", col=col, 
       main=paste("MA plot:", comparison, sep="\n"), ylim=lim, 
       xlab="Average expression", ylab="logFC")
  abline(h=0,col="darkgrey", lty=3)
}


WriteToHTML <- function (top, html.output, condition1, condition2) {
  c1.samples <- which(condition == condition1)
  c2.samples <- which(condition == condition2)
  top.de.genes <- as.character(top$ID)
  top.de.counts <- cbind(counts[top.de.genes,c1.samples],
                         counts[top.de.genes,c2.samples])
  top.de.norm.counts <- formatC(nc[top.de.genes,c(c1.samples,c2.samples)],
                               digits=2, format="f")

  # Combine raw and normalised counts in one cell
  html.count <- c()
  for (i in 1:nrow(top)) {
    html.count <- rbind(html.count, paste(top.de.counts[i,]," (", 
                        top.de.norm.counts[i,], ")", sep=""))
  }
  colnames(html.count) <- c(colnames(counts)[c1.samples],
                            colnames(counts)[c2.samples])

  # Add hyperlink to Ensembl ID column
  html.table <- cbind(top, html.count)
  html.table$ID <- sprintf("<a href=\"http://www.ensembl.org/id/%s\">%s</a>",
                           html.table$ID, html.table$ID)
  html.object <- xtable(html.table)

  # Add CSS to html
  css <- paste("<style type=\"text/css\">",
               "table {",
               "        border-width: 2px;",
               "        border-spacing: 2px;",
               "        border-style: solid;",
               "        border-color: gray;",
               "        border-collapse: collapse;",
               "}",
               "table td {",
               "        border-width: 1px;",
               "        padding: 4px;",
               "}",
               "</style>",
               sep="\n")
  tmp.filename <- "html.tmp"
  cat(css,file=tmp.filename)

  print.xtable(html.object, type="html", file=tmp.filename, 
               include.rownames=FALSE, append=TRUE)

  # Fix '>' and '<' in output
  sed.command <- sprintf("sed 's/&lt[ ;]/</g' %s | sed 's/&gt[ ;]/>/g' > %s; rm %s", tmp.filename, html.output, tmp.filename)
  system(sed.command)
}


TopHeatmap <- function(top, norm.counts, condition1, condition2) {
  c1.samples <- which(condition == condition1)
  c2.samples <- which(condition == condition2)
  # Select the top 20 genes by p-value
  if (annotate) {
    selected <- as.character(top$ID[1:20])
    gene.names <- paste(annotations[selected,1]," (",selected,")",sep="")
  } else {
    selected <- as.character(top$ID[1:20])
    gene.names <- rownames(nc[selected,])
  }
  log.top <- norm.counts[selected, c(c1.samples,c2.samples)]
  condition.labels <- c(rep(condition1, length(c1.samples)), 
                        rep(condition2, length(c2.samples)))
  sample.name.labels <- paste(condition.labels, colnames(log.top))
  covariate.labels <- condition.labels
  if (n.covariates > 0) {
    for (i in 1:n.covariates) {
     covariate.str <- sprintf("covariate.labels <- paste(covariate.labels, c%d[c(c1.samples, c2.samples)])", i)
     eval(parse(text=covariate.str))
    }
  }
  
  # Heatmap with covariate labels
  colnames(log.top) <- covariate.labels
  rownames(log.top) <- LETTERS[1:20]
  data.m <- melt(log.top)
  data.m <- ddply(data.m, .(X1), transform, rescaled = rescale(value))
  
  p <- ggplot(data.m, aes(X1, X2)) + 
       geom_tile(aes(fill = rescaled), colour="white") + 
       labs(x = "Top ranking genes", y = "Samples", 
            title="Heatmap of top 20 differentially expressed genes")
  g <- tableGrob(cbind(rownames(log.top), "Top ranking genes"=gene.names), 
                 gpar.coretext = gpar(fontsize=8),
                 gpar.coltext=gpar(fontsize=8),
                 gpar.rowtext=gpar(fontsize=8))

  grid.arrange(p + scale_fill_gradient(low = "white", high = "steelblue"),
               g, widths=c(6,2), nrow=1)

  # Heatmap with sample names
  colnames(log.top) <- sample.name.labels
  rownames(log.top) <- LETTERS[1:20]
  data.m <- melt(log.top)
  data.m <- ddply(data.m, .(X1), transform, rescaled = rescale(value))

  p <- ggplot(data.m, aes(X1, X2)) + 
       geom_tile(aes(fill = rescaled), colour="white") + 
       labs(x = "Top ranking genes", y = "Samples", 
            title="Heatmap of top 20 differentially expressed genes")
  g <- tableGrob(cbind(rownames(log.top), "Top ranking genes"=gene.names), 
                 gpar.coretext = gpar(fontsize=8),
                 gpar.coltext=gpar(fontsize=8), 
                 gpar.rowtext=gpar(fontsize=8))
  
  grid.arrange(p + scale_fill_gradient(low = "white", high = "steelblue"),
               g, widths=c(6,2), nrow=1)
}


MultipleMDSPlots <- function(data, mds.output, max.dimensions=8,
                             selection="pairwise", numPositions=1000) {
  colours <- rainbow(n=length(unique(condition)),end=0.7)
  col <- colours[as.numeric(condition)]
  pdf(file=mds.output)
  for (i in 1:min(ncol(counts)-2, max.dimensions-1)) {
    j <- i + 1
    plotMDS(data, top=numPositions, labels=colnames(data), col=col, 
            gene.selection=selection, dim.plot=c(i,j), cex=0.6, 
            main=sprintf("MDS plot using top %d genes", numPositions))
    legend("bottomleft", legend=levels(condition), col=colours, pch=19)
  }
  dev.off()
}

TopHistograms <- function(top, norm.counts, condition1, condition2,
                          n=20) {
  c1.samples <- which(condition == condition1)
  c2.samples <- which(condition == condition2)
  c.samples <- c(c1.samples, c2.samples)
  c.labels <- c(rep(condition1, length(c1.samples)),
                rep(condition2, length(c2.samples)))
  sample.order <- factor(colnames(norm.counts)[c.samples], 
                         levels=colnames(norm.counts)[c.samples])
  # Select the top 20 genes by p-value
  n <- ceiling(n/4) * 4
  selected <- as.character(top$ID[1:n])
  if (annotate) {
    gene.names <- paste(annotations[selected,1]," (",selected,")",sep="")
  } else {
    gene.names <- rownames(nc[selected,])
  }
  p <- list()
  for (i in 1:length(selected)) {
    j <- ceiling(i/4)
    k <- ifelse(test=(i %% 4 == 0), yes=4, no=i %% 4)
    df <- data.frame(sm=sample.order, condition=c.labels,
                     nc=norm.counts[selected[i], c.samples])
    if (i %% 4 == 1) {
      p[[j]] <- list()
    }
    p[[j]][[k]] <- ggplot(data=df, aes(x=sm, y=nc)) + 
                   geom_bar(aes(fill=condition), stat="identity") +
                   theme(axis.text.x=element_text(angle=90),
                         axis.text=element_text(size=8)) +
                   labs(x="", y="Normalised expression", title=gene.names[i])
    if (i %% 4 == 0) {
      do.call("grid.arrange", p[[j]])
    }
  }
}

PlotPCA <- function(data, pca.output, numPositions=1000) {
  if (ncol(counts) > 4) {
    # Use numPositions with the most variance  
    vars <- apply(data,1,var)
    positions <- names(sort(vars, decreasing=TRUE)[1:numPositions])
    dat <- t(data[positions,])
    pca <- prcomp(dat)
    pal <- rainbow(n=length(unique(condition)),end=0.7)
    col <- pal[as.numeric(condition)]
    pdf(file=pca.output)
    par(xpd=FALSE)
    pairs(pca$x[,1:4],col=col, pch=19, oma=c(4,4,6,20), cex=0.5)
    par(xpd=TRUE)
    legend(0.85, 0.55, legend=unique(condition), col=unique(col), pch=19)
    dev.off()
 } else {
    print("Not enough samples for PCA plots.")
  }
}


TopVarianceHeatmap <- function(data, heatmap.output, numPositions=1000) {
  # Use numPositions with the most variance
  vars <- apply(data,1,var)
  positions <- names(sort(vars, decreasing=TRUE)[1:numPositions])
  dat <- data[positions,]
  pdf(file=heatmap.output)
  par(cex.main=0.5)
  heatmap(dat, labCol=sample.names, 
          margins=c(9,3), cexCol=0.5, cexRow=0.5, cex.main=0.5, labRow=F,
          main=sprintf("Heatmap of %d genes with the most variance", 
          numPositions))
  dev.off()
}


######################################################################
## Load data

load(count.rdata)

# Check if annotations available
if (length(annotations) == 0) {
  annotate <- FALSE
} else {
  annotate <- TRUE
}

# Make a copy of the csv files in the output directory
system(paste("cp",sample.csv.filename,"samples.csv"))
system(paste("cp",comparison.csv.filename,"comparisons.csv"))


#####################################################################
## Create DGE list and filter genes

y <- DGEList(counts=counts,group=condition)
y$samples

# Filter genes. Only include genes where expression in cpm (counts per million) #   is greater than 1 in at least two samples
isexpr <- rowSums(cpm(y)>1) >= 2
y <- y[isexpr,]
y$samples$lib.size <- colSums(y$counts)

# How many fail and how many pass the filter
table(isexpr)
dim(y)

# Calculate normalisation factors of count matrix
y <- calcNormFactors(y, method="TMM")
y$samples

log.counts <- log(y$counts+1)
nc <- cpm(y, normalized.lib.sizes=T)
log.nc <- log(nc + 1)


#####################################################################
## Boxplots of gene counts

pdf(file="plots.pdf")
par(mar=c(9.1,4.1,4.1,2.1))

# Boxplots of log(raw) and log(TMM normalised) gene counts
# NB. Outliers are hidden
boxplot(log.counts, outline=F, ylab="log raw gene count", 
        main="Log raw counts", las=2, cex.axis=0.5)

boxplot(log.nc, outline=F, ylab="log TMM normalised gene count", 
        main="Log TMM Normalised", las=2, cex.axis=0.5)


#####################################################################
## RLE plots

gene.medians <- apply(log.counts,1,median)
boxplot(log.counts - gene.medians, outline=F, ylab="RLE",
        main="RLE plot of raw counts", las=2, cex.axis=0.5)
abline(h=0, lty=3)

gene.medians <- apply(log.nc,1,median)
boxplot(log.nc - gene.medians, outline=F, ylab="RLE",
        main="RLE plot of TMM-normalised counts", las=2, cex.axis=0.5)
abline(h=0, lty=3)

par(mar=c(5.1,4.1,4.1,2.1))


#####################################################################
## Create design matrix

if (n.covariates > 0) {
  design.str <- "design <- model.matrix(~0+condition"
  for (i in 1:n.covariates) {
     covariate.str <- sprintf("c%d <- factor(covariates[,%d])", i, i)
     eval(parse(text=covariate.str))
     design.str <- paste(design.str, sprintf("+c%d", i), sep="")
  }
  design.str <- paste(design.str, ")", sep="")
  print(design.str)
  eval(parse(text=design.str))
} else {
  design <- model.matrix(~0+condition)
}
colnames(design)[1:length(levels(condition))] <- levels(condition)
design


#####################################################################
## Create comparison matrix

cm.str <- "cm <- makeContrasts("
for (i in 1:n.comparisons) {
  condition1 <- comparisons[i,1]
  condition2 <- comparisons[i,2]
  cm.str <- paste(cm.str, condition1, "_vs_", condition2, "=", condition1, 
                  "-", condition2, ",", sep="")
}
cm.str <- paste(cm.str, "levels=design)", sep="")
print(cm.str)
eval(parse(text=cm.str))
cm


#####################################################################
## Run Voom or edgeR

if (voom) {

  ###################################################################
  ## Run Voom and plot plots
  
  print("Running Voom analysis")
  v <- voom(y, design, plot=TRUE)
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  dev.off()
  
  print(summary(decideTests(fit2)))
  
  for (comparison in colnames(cm)) {
    top.output <- paste("top_genes_",comparison,".txt",sep="")
    all.output <- paste("all_genes_",comparison,".txt",sep="")
    html.output <- paste("top_genes_",comparison,".html",sep="")
    pdf(paste("plots_",comparison,".pdf",sep=""), height=7, width=10)

    condition1 <- names(which(cm[,comparison] != 0))[1]
    condition2 <- names(which(cm[,comparison] != 0))[2]
    
    # MA plot
    VoomPlotMA(fit2, comparison)
    
    # P-value distribution
    hist(topTable(fit2, coef=comparison, adjust="BH", sort="p", n=Inf)$P.Value, 
         breaks=seq(0,1,0.02), main=paste("P-Value Distribution\n", 
         comparison, sep=""), xlab="P-value")

    # Print out top 10 DE genes to stdout
    print(topTable(fit2, coef=comparison, adjust="BH",sort="p", n=10), 
          row.names=FALSE)

    # Write top genes to file with annotations
    top <- format(topTable(fit2, coef=comparison, n=n.top, 
                  adjust="BH", sort="p"), digits=4, trim=TRUE)  
    if (is.null(top$ID)) {
      top <- cbind(ID=rownames(top), top)
      top$ID <- as.character(top$ID)
    }
    if (annotate) {
      top <- cbind(top, annotations[top$ID,])
      top$ID <- as.character(top$ID)
    }
    write.table(top, file=top.output, quote=FALSE, row.names=FALSE, sep="\t")

    # Write all genes to text file (no annotations)
    all <- format(topTable(fit2, coef=comparison, n=Inf, 
                  adjust="BH", sort="p"), digits=4, trim=TRUE)
    if (is.null(all$ID)) {
      all <- cbind(ID=rownames(all), all)
      all$ID <- as.character(all$ID)
    }
    write.table(all, file=all.output, 
                quote=FALSE, row.names=FALSE, sep="\t")
    
    # Write top genes to html file
    WriteToHTML(top, html.output, condition1, condition2)
    
    # Plot heatmap
    TopHeatmap(top, v$E, condition1, condition2)
    
    # Plot histograms
    TopHistograms(top, log.nc, condition1, condition2)
    dev.off()
  }
  
  # Plot MDS for up to 8 dimensions
  MultipleMDSPlots(v$E, "MDS.pdf", 8, "pairwise", numPositions=1000)
  
  # Plot Pair PCA plots
  PlotPCA(v$E, "PCA.pdf", numPositions=1000)

  # Clustered heatmap  
  TopVarianceHeatmap(v$E, "heatmap.pdf", numPositions=1000)
    
} else {
  
  ###################################################################
  ## Run edgeR and plot plots
  
  print("Running edgeR analysis")
  # estimate dispersion across genes
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  
  # visualise dispersion
  plotBCV(y)
  dev.off()
  
  # fit linear model
  fit <- glmFit(y, design)
  
  for (comparison in colnames(cm)) {
    top.output <- paste("top_genes_",comparison,".txt",sep="")
    all.output <- paste("all_genes_",comparison,".txt",sep="")
    html.output <- paste("top_genes_",comparison,".html",sep="")
    pdf(paste("plots_",comparison,".pdf",sep=""), height=7, width=10)
    
    condition1 <- names(which(cm[,comparison] != 0))[1]
    condition2 <- names(which(cm[,comparison] != 0))[2]
    
    # likelihood ratio test for condition
    lrt <- glmLRT(fit, contrast=cm[,comparison])
  
    # MA plot
    o <- order(lrt$table$PValue)
    cpm(y)[o[1:10],]
  
    # How many DE genes at 5% FDR
    print(summary(de <- decideTestsDGE(lrt)))
  
    # MA plot
    detags <- rownames(y)[as.logical(de)]
    plotSmear(lrt, de.tags=detags)
    abline(h=c(-1, 1), col="blue")
  
    # plot P-value distribution
    hist(topTags(lrt, n=Inf, sort.by="p")$table$PValue, 
         breaks=seq(0,1,0.025), main="P-Value Distribution")
  
    # Print out top 10 DE genes to stdout
    print(topTags(lrt, adjust.method="BH", sort="p", n=10), row.names=FALSE)

    # Write top genes to file with annotations
    top <- format(topTags(lrt, n=n.top, adjust.method="BH", sort.by="p")$table, 
                  digits=4, trim=TRUE)
    if (is.null(top$ID)) {
      top <- cbind(ID=rownames(top), top)
    }
    if (annotate) {
      top <- cbind(top, annotations[as.character(top$ID),])
    }
    write.table(top, file=top.output, quote=FALSE, row.names=FALSE, sep="\t")

    # Write all genes to text file (no annotations)
    all <- format(topTags(lrt, n=Inf, adjust.method="BH", sort.by="p")$table, 
                  digits=4, trim=TRUE)
    all <- cbind(ID=rownames(all), all)
    write.table(all, file=all.output,
                quote=FALSE, row.names=FALSE, sep="\t")
    
    # Write top genes to html file
    WriteToHTML(top, html.output, condition1, condition2)
    
    # Plot heatmap
    TopHeatmap(top, log.nc, condition1, condition2)
    
    # Plot histograms
    TopHistograms(top, log.nc, condition1, condition2)
    dev.off()

  }
  
  # Plot MDS for up to 8 dimensions
  MultipleMDSPlots(log.nc, "MDS.pdf", 8, "pairwise", numPositions=1000)
  
  # Plot Pair PCA plots
  PlotPCA(log.nc, "PCA.pdf", numPositions=1000)

  # Clustered heatmap  
  TopVarianceHeatmap(log.nc, "heatmap.pdf", numPositions=1000)

}


######################################################################
## Save image

# # uncomment to save image
# save.image("image.RData")


proc.time()
sessionInfo()