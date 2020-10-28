library(parallel)
library(Matrix)
library(edgeR)

loadRData <- TRUE
compareDiseases <- FALSE
findMarkers <- TRUE
gene_freq_cutoff <- 0.1  # cutoff to decide how many genes are tested
nCore <- 15  # number of cores for mclapply


source("functions.R")

if (loadRData == FALSE){  # requires > 100GB memory.
    timeStamp("Data loading...")

    counts <- readMM("counts_QCed_1000gene10MT.mtx.gz")
    anno <- read.csv("annotation_QCed_1000gene10MT.csv", stringsAsFactors=FALSE)
    ecm <- read.csv("ECM_matrisome_hs_masterlist.csv", stringsAsFactors=FALSE)
    ecm <- ecm[, c("Division", "Category", "Gene.Symbol")]
    genes <- read.csv("geneNames_nonZeroRemoved.csv", stringsAsFactors=FALSE)
    cells <- read.csv("cellNames_QCed_1000gene10MT.csv", stringsAsFactors=FALSE)
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
    rm(notCOPD, notMultiplet, notPNEC, expressedGenes)

    timeStamp("Log2-cpm transform the data.")
    # since the data are huge, transform them by batch and concatenate in the end.
    batch_size <- round(dim(counts)[2] / 15)
    col_indices <- mclapply(
        mc.cores=nCore,
        1:15,
        function(i){
            range.begin <- 1 + batch_size * (i - 1)
            range.end <- batch_size * i
            if (i != 15){
                range.sub <- range.begin:range.end
            } else {range.sub <- range.begin:dim(counts)[2]}

            return(range.sub)
        }
    )
    list.counts.normalized <- mclapply(
        mc.cores=20,
        1:15,
        function(i){
            timeStamp(paste("Normalizing batch", i, "/ 15."))
            out <- cpm(counts[, col_indices[[i]]])
            out <- log2(out + 1)
            return(signif(out, 4))
        }
    )
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
        mc.cores=nCore,  # use max core - 1.
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




if (compareDiseases){
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

    timeStamp("T-test comparing disease vs normal done.")
}

if (findMarkers){
    # divide macrophages by disease states for marker gene identification if needed.
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


    timeStamp("Subset the count matrix by ECM genes.")
    counts.logcpm <- counts.logcpm[rownames(counts.logcpm) %in% ecm$Gene.Symbol, ]
    

    timeStamp("Find marker ECM genes in each cell type.")
    t.markers.list <- sapply(
        simplify=FALSE,
        unique(anno$Manuscript_Identity),
        function(cat){
            timeStamp(paste("T-test run for", cat, "to find cell type markers."))
            print(paste(which(unique(anno$Manuscript_Identity) == cat), "out of", length(unique(anno$Manuscript_Identity)), "categories."))

            t <- t_test.mc(counts.logcpm, anno$Manuscript_Identity == cat, cat)
            print(head(t))
            return(t)
        }
    )
    #timeStamp("T-test run for Macrophage to find cell type markers.")
    #t.markers.list[["Macrophage"]] <- t <- t_test(counts.logcpm, anno$Manuscript_Identity == "Macrophage", "Macrophage")

    # save the results.
    mclapply(
        mc.cores=nCore,
        names(t.markers.list),
        function(n){
            write.csv(t.markers.list[[n]], paste0("DE/tTest_markers_", n, ".csv"))
            return(NULL)
        }
    )
    # also save a combined result.
    write.csv(do.call(rbind, t.markers.list), "DE/tTest_markers_combined.csv")

    timeStamp("T-test to find cell type markers done.")
}

timeStamp("Done!")