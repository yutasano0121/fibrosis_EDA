library(Matrix)
library(ggplot2)

savePlots <- FALSE

# cutoff threshold for QC
cutoff.geneNum <- 1000
cutoff.MTratio <- 10


read_tsv <- function(f, HEADER=TRUE){
    out <- read.csv(f, sep="\t", header=HEADER)
    return(out)
}

# count matrix
counts <- readMM(file="GSE136831_RawCounts_Sparse.mtx.gz")
# cell ID
cells <- read_tsv("GSE136831_AllCells.cellBarcodes.txt.gz", HEADER=FALSE)
# gene ID
genes <- read_tsv("GSE136831_AllCells.GeneIDs.txt.gz")
# metadata
anno <- read_tsv("GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")

rownames(counts) <- genes$HGNC_EnsemblAlt_GeneID
colnames(counts) <- cells$V1
rm(genes, cells)


# number of detected genes
geneNum <- colSums(counts != 0)
anno$gene.number <- geneNum
print(paste("Mean gene count:", mean(geneNum)))

# MT gene ratio (%)
MT <- grep("MT-", rownames(counts), perl=TRUE, value=TRUE)
MTratio <- colSums(counts[MT,]) / colSums(counts) * 100
anno$MTratio <- MTratio
print(paste("Mean MT-gene ratio (%):", mean(MTratio)))


# remove cells with zero counts for endogenous genes.
notZeroCounts <- geneNum != 0
counts <- counts[, notZeroCounts]
anno <- anno[notZeroCounts, ]

# plot QC criteria
p.qc <- ggplot(anno, aes(x=MTratio, y=gene.number, color=Disease_Identity)) + geom_point()

if (savePlots){
    png("qc.png", height=5, width=6, unit="in", res=300)
    p.qc
    dev.off()
}

# set QC thresholds
qc <- (geneNum >= cutoff.geneNum) & (MTratio <= cutoff.MTratio)
print(paste0("Cutoff thresholds for QC: gene counts...", cutoff.geneNum, ", MT gene ratio (%)...", cutoff.MTratio))
print(paste("Number of cells which survived QC:", sum(qc), "out of", length(qc), "cells."))

# subset the data by the QC thresholds
counts <- counts[, qc]
anno <- anno[qc, ]

# remove unexpressed genes
counts <- counts[rowSums(counts) != 0, ]

# save QC'ed data
# writeMM doesn't save row/column names, so they need to be saved separately.
writeMM(counts, "counts_QCed_1000gene10MT.mtx")
write.csv(anno, "annotation_QCed_1000gene10MT.csv", row.names=FALSE)
write.csv(rownames(counts), "geneNames_nonZeroRemoved.csv")
write.csv(colnames(counts), "cellNames_QCed_1000gene10MT.csv")
print("QC'ed data saved.")

# randomly sample 10% of the data
set.seed(0)
downsample <- sample(x=c(1:dim(counts)[2]), size=round(dim(counts)[2]/10))
counts <- counts[, downsample]
anno <- anno[downsample, ]

# save the reduced data
writeMM(counts, "counts_QCed_1000gene10MT_reduced.mtx")
write.csv(colnames(counts), "cellNames_QCed_1000gene10MT_reduced.csv")
write.csv(anno, "annotation_QCed_1000gene10MT_reduced.csv", row.names=FALSE)
print("Downsampled data saved.")

print("All done!")

