library(parallel)
library(Matrix)
library(edgeR)

detectCores()

loadRData <- TRUE

source("functions.R")

if (loadRData == FALSE){  # requires > 100GB memory.
    timeStamp("Data loading...")

    counts <- readMM("counts_QCed_1000gene10MT.mtx.gz")
    anno <- read.csv("annotation_QCed_1000gene10MT.csv", stringsAsFactors=FALSE)
    genes <- read.csv("geneNames_nonZeroRemoved.csv")
    cells <- read.csv("cellNames_QCed_1000gene10MT.csv")
    ecm <- read.csv("ECM_matrisome_hs_masterlist.csv", stringsAsFactors=FALSE)
    timeStamp("Data loaded.")

    rownames(counts) <- genes$x
    colnames(counts) <- cells$x
    rm(genes, cells)

    timeStamp("Remove COPD, multiplets, PNECs, and ionocytes.")
    # remove PNEC and ionocytes, since it has only 22 and 19 observation.
    notCOPD <- anno$Disease_Identity != "COPD"
    notMultiplet <- anno$CellType_Category != "Multiplet"
    notPNEC <- anno$Manuscript_Identity != "PNEC"
    notIonocyte <- anno$Manuscript_Identity != "Ionocyte"
    counts <- counts[, (notCOPD & notMultiplet) & (notPNEC & notIonocyte)]
    anno <- anno[(notCOPD & notMultiplet) & (notPNEC & notIonocyte), ]

    timeStamp("Remove genes not expressed.")
    expressedGenes <- rowSums(counts) != 0
    counts <- counts[expressedGenes, ]
    ecm <- ecm[ecm$Gene.Symbol %in% rownames(counts), ]
    rm(notCOPD, notMultiplet, notPNEC, expressedGenes)

    timeStamp("Log2-cpm transform the data.")
    # since the data are huge, transform them by batch and concatenate in the end.
    batch <- round(dim(counts)[2] / 5)
    range.begin <- 1
    range.end <- batch
    col_indices <- list()
    for (i in 1:5){
        if (i != 5){
            range.sub <- range.begin:range.end
        } else {range.sub <- range.begin:dim(counts)[2]}

        col_indices[[i]] <- range.sub

        range.begin <- range.end + 1
        range.end <- range.end + batch
    }
    list.counts.normalized = list()
    for (i in 1:5){
        timeStamp(paste0("Normalizing batch ", i, "/5."))
        list.counts.normalized[[i]] <- cpm(counts[, col_indices[[i]]])
    }
    rm(counts)

    # combine them into a single data frame.
    counts.logcpm <- do.call(cbind, list.counts.normalized)
    rm(list.counts.normalized)
    timeStamp("Counts log2-cpm transformed.")

    save(counts.logcpm, anno, ecm, file="RData/QCed_log2cpm_allCells_allGenes.RData")

    timeStamp("Data loaded and normalized. RData saved.")
} else {
    timeStamp("RData loading...")
    load("RData/QCed_log2cpm_allCells_allGenes.RData")
    timeStamp("RData loaded.")
}
rm(ecm)

gene_freq_cutoff <- 0.1  # cutoff to decide how many genes are tested


timeStamp("Print distribution of IPF/Control in each category.")
cellTypeDistribution <- lapply(
    unique(anno$Manuscript_Identity),
    function(cat){
        num_IPF <- sum((anno$Manuscript_Identity == cat) & (anno$Disease_Identity == "IPF"))
        num_total <- sum(anno$Manuscript_Identity == cat)
        print(paste(
            cat, " cells have", num_total, "cells in total, in which", 
            num_IPF, "cells,", num_IPF / num_total * 100, "% are IPF."
        ))
        return(c(cellType=cat, total=num_total, IPF=num_IPF, IPFratio=num_IPF/num_total))
    }
)
# write.csv(do.call(rbind, cellTypeDistribution), "cellTypeDistribution.csv")


t_test.mc <- function(DF, CATEGORY, CATEGORY_NAME=NULL){  # DF = genes x cells
    # remove genes not expressed by the "TRUE" group
    genes.expressed <- find_expressedGenes(DF[, CATEGORY], gene_freq_cutoff)
    print(paste("T-test is run for", sum(genes.expressed), "genes."))

    t.result <- mclapply(
        mc.cores=15,  # use max core - 1.
        rownames(DF)[genes.expressed],
        function(gene){
            print(paste(
                which(rownames(DF)[genes.expressed]==gene), "/", sum(genes.expressed), "genes", 
                which(unique(anno$Manuscript_Identity) == CATEGORY_NAME), "/", length(unique(anno$Manuscript_Identity)), "categories"
            ))
            t <- t.test(DF[gene, ]~CATEGORY)
            if (is.null(CATEGORY_NAME)){
                out <- c(
                    gene=gene,
                    logFC=t$estimate[[2]] - t$estimate[[1]],  # TRUE - FALSE
                    p.value=t$p.value
                )    
            } else {
                out <- c(
                    cellType=CATEGORY_NAME,
                    gene=gene,
                    logFC=t$estimate[[2]] - t$estimate[[1]],  # TRUE - FALSE
                    p.value=t$p.value
                )
            }
            return(out)
        }
    )

    # make it a data frame
    t.result <- as.data.frame(do.call(rbind, t.result), stringsAsFactors=FALSE)
    # make logFC and p.values numeric
    t.result$logFC <- as.numeric(t.result$logFC)
    t.result$p.value <- as.numeric(t.result$p.value)
    # calculate FDR
    t.result$fdr <- p.adjust(t.result$p.value, "BH")
    # remove non-significant or downregulated genes
    t.result <- t.result[(t.result$fdr < 0.05) & (t.result$logFC > 0), ]
    # sort by p-values
    t.result <- t.result[order(t.result$fdr), ]

    print(paste(dim(t.result)[1], "DEGs found."))

    return(t.result)
}





timeStamp("Compare each cell subtypes for IPF vs. Control.")
categories <- unique(anno$Manuscript_Identity)
cat_exclude <- c("Aberrant_Basaloid")
categories <- categories[!categories %in% cat_exclude]
t.IPFvsCtrl.list <- sapply(
    simplify=FALSE,
    categories,
    function(cat){
        timeStamp(paste("T-test run for", cat, "comparing IPF vs. Control."))
        print(paste(which(categories == cat), "out of", length(categories), "categories."))

        # subset by category.
        cat_to_test <- anno$Manuscript_Identity == cat
        # compare IPF vs. Control.
        t <- t_test.mc(counts.logcpm[, cat_to_test], anno[cat_to_test, ]$Disease_Identity == "IPF", cat)
        print(head(t))
        return(t)
    }
)
# compare Abberant Basaloids with other epithelial cells.
t.IPFvsCtrl.list[["Aberrant_Basaloid"]] <- t_test.mc(
    counts.logcpm[, anno$CellType_Category=="Epithelial"], 
    anno[anno$CellType_Category=="Epithelial", ]$Manuscript_Identity == "Aberrant_Basaloid", 
    "Aberrant_Basaloid"
)

# save the results.
lapply(
    names(t.IPFvsCtrl.list),
    function(n){
        write.csv(t.IPFvsCtrl.list[[n]], paste0("DE/tTest_IPFvsCtrl_", n, ".csv"))
        return(NULL)
    }
)
# also save a combined result.
write.csv(do.call(rbind, t.IPFvsCtrl.list), "DE/tTest_IPFvsCtrl_combined.csv")

if (FALSE){
    # divide macrophages by disease states for marker gene identification.
    Mf_d <- (anno$Manuscript_Identity == "Macrophage") & (anno$Disease_Identity == "IPF")
    Mf_n <- (anno$Manuscript_Identity == "Macrophage") & (anno$Disease_Identity == "Control")
    MfA_d <- (anno$Manuscript_Identity == "Macrophage_Alveolar") & (anno$Disease_Identity == "IPF")
    MfA_n <- (anno$Manuscript_Identity == "Macrophage_Alveolar") & (anno$Disease_Identity == "Control")
    anno$category <- anno$Manuscript_Identity
    anno$category[Mf_d] <- "Macrophage_IPF"
    anno$category[Mf_n] <- "Macrophage_Control"
    anno$category[MfA_d] <- "Mf_Alveolar_IPF"
    anno$category[MfA_n] <- "Mf_Alveolar_Control"
    rm(Mf_d, Mf_n, MfA_d, MfA_n)


    timeStamp("Find marker genes in each category.")
    t.markers.list <- sapply(
        simplify=FALSE,
        unique(anno$category),
        function(cat){
            timeStamp(paste("T-test run for", cat, "to find cell type markers."))
            print(paste(which(unique(anno$category) == cat), "out of", length(unique(anno$category)), "categories."))

            # compare IPF vs. Control.
            t <- t_test(counts.logcpm, anno$category == cat, cat)
            print(head(t))
            return(t)
        }
    )
    timeStamp("T-test run for MAcrophage to find cell type markers.")
    t.markers.list[["Macrophage"]] <- t <- t_test(counts.logcpm, anno$Manuscript_Identity == "Macrophage", "Macrophage")

    # save the results.
    lapply(
        names(t.markers.list),
        function(n){
            write.csv(t.markers.list[[n]], paste0("DE/tTest_markers_", n, ".csv"))
            return(NULL)
        }
    )
    # also save a combined result.
    write.csv(do.call(rbind, t.markers.list), "DE/tTest_markers_combined.csv")
}

timeStamp("Done!")