library(Matrix)
library(edgeR)
library(Rtsne)
library(umap)
library(ggplot2)
library(viridis)

if (FALSE){
    cpm <- function(df){
        out <- apply(
            df, 2,
            function(x){
                return(x / sum(x) * 1E6)
            }
        )
        return(out)
    }
}

source("functions.R")
timeStamp("Data loading...")
counts <- readMM("counts_QCed_1000gene10MT_reduced.mtx.gz")
anno <- read.csv("annotation_QCed_1000gene10MT_reduced.csv")
genes <- read.csv("geneNames_nonZeroRemoved.csv")
cells <- read.csv("cellNames_QCed_1000gene10MT_reduced.csv")
ecm <- read.csv("ECM_matrisome_hs_masterlist.csv", stringsAsFactors=FALSE)
timeStamp("Data loaded.")

rownames(counts) <- genes$x
colnames(counts) <- cells$x

# remove genes not expressed
expressedGenes <- rowSums(counts) != 0
counts <- counts[expressedGenes, ]

# log2-cpm transform the data.
counts.logcpm <- edgeR::cpm(counts, log=TRUE, prior.count=1)
timeStamp("Counts log2-cpm transformed.")

# subset the data with ECM-associated genes.
counts.ecm <- counts.logcpm[rownames(counts.logcpm) %in% ecm$Gene.Symbol, ]


# mean-center the data.
# returns cells x genes
counts.ecm <- apply(
    counts.ecm, 1, function(row){
        return(row - mean(row))
    }
)

timeStamp("Start t-SNE.")
set.seed(0)
tsne <- Rtsne(counts.ecm)  
timeStamp("t-SNE done.")

p.celltype <- plot.tsne(DF=as.data.frame(tsne$Y), ANNO=anno)
p.disease <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Disease_Identity', COLOR.LAB='Diseases'
)
p.subject <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Subject_Identity', COLOR.LAB='Subjects'
)

pngStore(p.celltype, "tsne_reduced_cellType.png", W=10, H=10)
pngStore(p.disease, "tsne_reduced_disease.png", W=10, H=10)
pngStore(p.subject, "tsne_reduced_subject.png", W=15, H=10)

timeStamp("Done!")
