#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("Usage: DESeq2.r <merged_gene_counts.txt> <design CSV> <FDR cutoff> <log2 fold change cutoff> <heatmap column order> <comparisons CSV(optional)>", call.=FALSE)
}
merged_counts <- args[1]
design_csv <- args[2]
fdr_cutoff <- as.numeric(args[3])
lfc_cutoff <- as.numeric(args[4])
heatmap_group_order <- if(args[5]=="true") TRUE else FALSE

# Load / install packages
if (!require("DESeq2")) {
    install.packages("BiocManager")
    BiocManager::install("DESeq2")
    BiocManager::install("limma")
    library("DESeq2")
    library("limma")
}
if (!require("pheatmap")) {
    install.packages("pheatmap")
    library("pheatmap")
}
if (!require("ggplot2")) {
    install.packages("ggplot2")
    library("ggplot2")
}
if (!require("openxlsx")) {
    install.packages("openxlsx")
    library("openxlsx")
}

# Read the merged counts
countdata <- read.csv(merged_counts, header=TRUE, row.names=1, sep=",", check.names=FALSE, stringsAsFactors=FALSE)
# Separate gene names
gene_name <- data.frame(countdata$gene_name, row.names = row.names(countdata), check.names=FALSE, stringsAsFactors=FALSE)
# Read the sample grouping
grouping <- read.csv(design_csv, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# Only use samples that appear in both design CSV and merged counts
shared_snames <- intersect(grouping$sample, names(countdata))
grouping <- grouping[grouping$sample %in% shared_snames,]
countdata <- as.matrix(countdata[,grouping$sample, drop=FALSE])
# If only one group label present, set it to NA to avoid DESeq2 error
if (length(unique(grouping$group)) < 2) { grouping$group <- NA }
# If group label absent, use sample label
grouping$group <- ifelse(is.na(grouping$group), grouping$sample, grouping$group)
# Filter groups so that all groups have at least 2 replicates
no.replicates <- table(grouping$group)
groups <- names(no.replicates)[no.replicates > 1]

# Create DESeq object
# Create coldata
coldata <- data.frame(row.names=grouping$sample, group=grouping$group)
if (ncol(countdata) > 1) {
    # Load DESeq object
    dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group)
} else {
    # When there is only one sample, design has to be ~1
    dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
}
# Drop rows with all zeroes
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- estimateSizeFactors(dds)
# Export normalized counts
normalized_counts <- counts(dds,normalized=TRUE)
write.xlsx(data.frame("Gene ID"=row.names(normalized_counts),normalized_counts, check.names=FALSE), "normalized_counts.xlsx")
# Transform the data using rlog or log2(x+1) depending on whether there is a group with replicates
if (length(groups) > 0) {
    transformed <- rlog(dds, blind=FALSE)
} else {
    transformed <- normTransform(dds)
    coldata <- NA # No need to label figures with group labels if there are no replicates
}

# Make study-level figures, only happens if there are more than 2 samples
if (length(grouping$sample) > 2) {
    # MDS plot
    jpeg("DESeq2_sample_MDS_plot.jpg", width=8, height=8, unit="in", res=300)
    p <- limma::plotMDS(assay(transformed))
    dev.off()
    # Also export plot data for MultiQC
    # mdsdata <- data.frame(p$cmdscale.out)
    # names(mdsdata) <- c("Dimension 1", "Dimension 2")
    # mdsdata$group <- ifelse(is.na(coldata), NA, coldata$group)
    # write.table(mdsdata, "DESeq2_sample_MDS_plot.tsv", sep="\t", quote=FALSE)
    # Extract MDS coordinates
    mdsdata <- data.frame(Dimension1 = p$x, Dimension2 = p$y)
    mdsdata$group <- ifelse(is.na(coldata), NA, coldata$group)
    write.table(mdsdata, "DESeq2_sample_MDS_plot.tsv", sep="\t", quote=FALSE)
    # Similarity matrix as heatmap
    correlation_matrix <- data.frame(cor(assay(transformed), method="pearson"))
    jpeg("DESeq2_sample_similarity_matrix.jpg", width=8, height=8, unit="in", res=300)
    p <- pheatmap(correlation_matrix, annotation_col=coldata, clustering_distance_cols="correlation", clustering_distance_rows="correlation")
    dev.off()
    # Also output the matrix to plot in MultiQC
    correlation_matrix <- correlation_matrix[p$tree_row$label[p$tree_row$order], p$tree_col$label[p$tree_col$order]]
    write.table(correlation_matrix, "DESeq2_sample_similarity_matrix.tsv", sep="\t", quote=FALSE)
    # Make a heatmap for top 100 genes with most variance
    topVarGenes <- head(order(rowVars(assay(transformed)), decreasing=TRUE), 100)
    topGeneMat <- assay(transformed)[topVarGenes, ]
    topGeneMat <- topGeneMat - rowMeans(topGeneMat)
    topGeneName <- gene_name[row.names(topGeneMat),]
    row.names(topGeneMat) <- ifelse(is.na(topGeneName), row.names(topGeneMat), topGeneName)
    if (heatmap_group_order) {
	topGeneMat <- topGeneMat[,row.names(coldata)[order(coldata$group)]]
    }
    jpeg("DESeq2_top_gene_heatmap.jpg", width=8, height=12, unit="in", res=300)
    pheatmap(topGeneMat, annotation_col=coldata, cluster_cols = !heatmap_group_order, fontsize_row=8)
    dev.off()
    # Also output the heatmap matrix to plot in MultiQC
    rowclust <- hclust(dist(topGeneMat))
    if(heatmap_group_order){
        topGeneMat <- topGeneMat[rowclust$order,] 
    }
    else {
    	colclust <- hclust(dist(t(topGeneMat)))
        topGeneMat <- topGeneMat[rowclust$order, colclust$order]
    }
    write.table(topGeneMat, "DESeq2_top_gene_heatmap.tsv", sep="\t", quote=FALSE)
}

# Carry out DGE analysis if there are at least 2 groups with replicates
if (length(groups) > 1) {
    # Remove groups with no replicates and carry out DEG analysis
    dds <- dds[,colData(dds)$group %in% groups]
    dds$group <- factor(dds$group)
    dds <- DESeq(dds)
    # Remove groups with no replicates and calculate geomean of transformed counts per group
    transformed <- transformed[,colData(transformed)$group %in% groups]
    transformed$group <- factor(transformed$group)
    transMeanPerLvl <- sapply(levels(transformed$group), function(lvl) rowMeans(assay(transformed)[,transformed$group==lvl]))
    transMeanPerLvl <- as.data.frame(2**transMeanPerLvl)
    # Read the group comparison file
    if (length(args) > 5) {
        comparisons_file <- args[6]
        comparisons <- read.csv(comparisons_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        comparisons <- t(comparisons)
        rownames(comparisons) <- c()
    } else {
        comparisons <- combn(groups, 2)
    }
    # Extract results of each pairs of groups
    for (i in 1:ncol(comparisons)) {
        cond1 <- comparisons[1,i]
        cond2 <- comparisons[2,i]
        if (!(cond1 %in% groups) || !(cond2 %in% groups)) {
            message(paste(cond1, "or", cond2, "did not have replicates. Skipping comparison..."))
            next
        }
        res <- results(dds, contrast=c("group",cond1,cond2), alpha=fdr_cutoff, lfcThreshold=lfc_cutoff)
        # Add transformed counts for each group
        res[,cond1] <- transMeanPerLvl[,cond1]
        res[,cond2] <- transMeanPerLvl[,cond2]
        # Shrink the fold change
        reslfc <- lfcShrink(dds, contrast=c("group",cond1,cond2), type="ashr")
        res$lfcShrink <- reslfc$log2FoldChange
        # Add the gene name back in
        res <- merge(gene_name, data.frame(res), by=0, all.y=TRUE)
        # When gene name is missing, default that to gene ID
        res[,2] <- ifelse(is.na(res[,2]), as.character(res[,1]), as.character(res[,2]))
        names(res)[1:2] <- c("Gene ID", "Gene Name")
        # Sort the results by FDR and pvalue
        res <- res[order(res$padj, res$pvalue), ]
        # Save results
        write.table(res, paste0(cond1,"_vs_",cond2,"_DESeq_results.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
        write.xlsx(res, paste0(cond1,"_vs_",cond2,"_DESeq_results.xlsx"))
        # Save MA plot
        jpeg(paste0(cond1,"_vs_",cond2,"_DESeq_MA_plot.jpg"), width=8, height=8, unit="in", res=300)
        plotMA(reslfc, alpha=fdr_cutoff)
        dev.off()
        # Plot scatter plot
        res <- as.data.frame(res[complete.cases(res),])
        plot_data <- res[c(cond1, cond2)]
        plot_data$DEG <- res$padj<fdr_cutoff
        jpeg(paste0(cond1,"_vs_",cond2,"_DESeq_scatterplot.jpg"), width=8, height=8, unit="in", res=300)
        print(ggplot(plot_data, aes_string(cond1, cond2, colour="DEG"))+geom_point()+scale_x_log10(limits=c(1,NA))+scale_y_log10(limits=c(1,NA)))
        dev.off()
    }
}