library(Matrix)
library(edgeR)
library(Rtsne)
library(umap)
library(ggplot2)
library(viridis)

saveRData <- TRUE
loadRData <- TRUE
removeStromal <- FALSE  # since they have too high expression of ECM genes...
removeTop <- TRUE  # remove top ECM expressors
ecmThreshold <- 0.975  # percentile threshold for ECM expression

# if edgeR is not used...
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

if (loadRData == FALSE){
    timeStamp("Data loading...")
    counts <- readMM("counts_QCed_1000gene10MT_reduced.mtx.gz")
    anno <- read.csv("annotation_QCed_1000gene10MT_reduced.csv")
    genes <- read.csv("geneNames_nonZeroRemoved.csv")
    cells <- read.csv("cellNames_QCed_1000gene10MT_reduced.csv")
    ecm <- read.csv("ECM_matrisome_hs_masterlist.csv", stringsAsFactors=FALSE)
    timeStamp("Data loaded.")

    rownames(counts) <- genes$x
    colnames(counts) <- cells$x

    # remove COPD cells.
    cells_notCOPD <- anno$Disease_Identity != "COPD"
    counts <- counts[, cells_notCOPD]
    anno <- anno[cells_notCOPD, ]

    # remove genes not expressed.
    expressedGenes <- rowSums(counts) != 0
    counts <- counts[expressedGenes, ]
    ecm <- ecm[ecm$Gene.Symbol %in% rownames(counts), ]

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


    if (saveRData==TRUE){
        save(counts.ecm, anno, tsne, file="RData/ecm_tsne_reduced.RData")
        timeStamp("RData saved.")
    }
} else {
    load("RData/ecm_tsne_reduced.RData")
    timeStamp("RData loaded.")
}


# make plots which do not require ECM information.
p.celltype <- plot.tsne(DF=as.data.frame(tsne$Y), ANNO=anno)
p.manuscript <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Manuscript_Identity', COLOR.LAB='Cell subtypes'
)
p.subclass <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Subclass_Cell_Identity', COLOR.LAB='Cell subtypes'
) + theme(legend.position="bottom")
p.disease <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Disease_Identity', COLOR.LAB='Diseases'
)
p.subject <- plot.tsne(
    DF=as.data.frame(tsne$Y), ANNO=anno,
    COLOR='Subject_Identity', COLOR.LAB='Subjects'
)


if (removeStromal == TRUE){
    print("Plot without stromal cells.")
    notEndo <- anno$CellType_Category != "Stromal"
    counts.ecm <- counts.ecm[notEndo, ]
    tsne$Y <- tsne$Y[notEndo, ]
    anno <- anno[notEndo, ]
} 


if (removeTop){
    print("Plot without top ECM expressors.")

    ecm.score <- sapply( # returns a named list.
        simplify=FALSE, 
        unique(ecm$Category),
        function(cat){
            genes.ecm.cat <- ecm$Gene.Symbol[ecm$Category==cat]
            score <- calcScore(t(counts.ecm), genes.ecm.cat)
        }
    )
    
    top.list <- lapply( # returna an unnamed list of cells to retain.
        names(ecm.score),
        function(cat){
            top <- ecm.score[[cat]] < quantile(ecm.score[[cat]], ecmThreshold)
            return(top) 
        }
    )

    # make a filter to remove top ECM expressors
    cells_to_retain <- rep(TRUE, dim(counts)[2])
    for (i in 1:length(top.list)){
        cells_to_retain <- cells_to_retain & top.list[[i]]
    }

    # remove top cells and re-scale the score.
    for (i in 1:length(ecm.score)){
        ecm.score[[i]] <- scale(ecm.score[[i]][cells_to_retain])
    }
    # convert it to a data frame.
    ecm.score <- as.data.frame(do.call(cbind, ecm.score))
    colnames(ecm.score) <- gsub(" ", "", unique(ecm$Category))

    counts.ecm <- counts.ecm[cells_to_retain, ]
    tsne$Y <- tsne$Y[cells_to_retain, ]
    anno <- anno[cells_to_retain, ]
} else {
    # calculate the score of ecm genes by categories.
    ecm.score <- sapply(
        simplify=FALSE, 
        unique(ecm$Category),
        function(cat){
            genes.ecm.cat <- ecm$Gene.Symbol[ecm$Category==cat]
            score <- calcScore(t(counts.ecm), genes.ecm.cat)
            return(scale(score))  # scale the score (scaled sum) again to obtain Z-scores!
        }
    )
    ecm.score <- as.data.frame(do.call(cbind, ecm.score))
    colnames(ecm.score) <- gsub(" ", "", unique(ecm$Category))
}


p.ecm.score <- sapply(
    simplify=FALSE,
    colnames(ecm.score),
    function(cat){
        p <- plot.tsne(
            DF=as.data.frame(tsne$Y), ANNO=ecm.score,
            COLOR=cat, COLOR.LAB=cat, 
            GENE=TRUE, GENE.DF=ecm.score
        )
    }
)
pngStore(p.celltype, "tsne_reduced_cellType.png", W=10, H=10)
pngStore(p.manuscript, "tsne_reduced_manuscriptIdent.png", W=12.5, H=10)
# pngStore(p.subclass, "tsne_reduced_subclass.png", W=10, H=15) 
pngStore(p.disease, "tsne_reduced_disease.png", W=10, H=10)
pngStore(p.subject, "tsne_reduced_subject.png", W=12.5, H=10)
if (removeStromal == TRUE){
    lapply(
        names(p.ecm.score),
        function(cat)
        pngStore(
            p.ecm.score[[cat]], 
            paste0("tsne_reduced_woStromal_", cat, ".png"), 
            W=10, H=10
        )
    )
} else if (removeTop == TRUE){
    lapply(
        names(p.ecm.score),
        function(cat)
        pngStore(
            p.ecm.score[[cat]], 
            paste0("tsne_reduced_woTop", gsub(" ", "", cat), ".png"), 
            W=10, H=10
        )
    )
} else {
    lapply(
        names(p.ecm.score),
        function(cat)
        pngStore(
            p.ecm.score[[cat]], 
            paste0("tsne_reduced_", gsub(" ", "", cat), ".png"), 
            W=10, H=10
        )
    )
}

timeStamp("Done!")
