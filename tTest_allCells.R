library(Matrix)
library(edgeR)

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

    timeStamp("Remove COPD, multiplets, and PNEC cells.")
    # remove PNEC, since it has only 22 observation.
    notCOPD <- anno$Disease_Identity != "COPD"
    notMultiplet <- anno$CellType_Category != "Multiplet"
    notPNEC <- anno$Manuscript_Identity != "PNEC"
    counts <- counts[, (notCOPD & notMultiplet) & notPNEC]
    anno <- anno[notCOPD & notMultiplet, ]

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

gene_freq_cutoff <- 0.05  # cutoff to decide how many genes are tested



# remove PNEC, since it has only 22 observation.
notPNEC <- anno$Manuscript_Identity != "PNEC"
counts.logcpm <- counts.logcpm[, notPNEC]
anno <- anno[notPNEC, ]


# divide macrophages by disease states
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


t_test <- function(DF, CATEGORY){
    # remove genes not expressed by each category
    genes.expressed1 <- find_expressedGenes(DF[, CATEGORY], gene_freq_cutoff)
    genes.expressed2 <- find_expressedGenes(DF[, !CATEGORY], gene_freq_cutoff)
    print(paste("T-test is run for", sum(genes.expressed1 & genes.expressed2), "genes."))
    DF <- DF[, genes.expressed1 & genes.expressed2]
    t.result <- sapply(
        rownames(DF),
        function(gene){
            t <- t.test(DF[gene, ]~CATEGORY)
            out <- c(
                p.value=t$p.value,
                logFC=t$estimate[[2]] - t$estimate[[1]]  # TRUE - FALSE
            )
            return(out)
        }
    )
    t.result <- as.data.frame(t(rbind(t.result)))
    # sort by p-values, and upregulated genes to the top
    t.result$fdr <- p.adjust(t.result$p.value, "BH")
    t.result <- t.result[order(t.result$fdr), ]
    t.result <- t.result[order(t.result$logFC < 0), ]
    # remove non-significant genes
    t.result <- t.result[t.result$fdr <= 0.05, ]

    return(t.result)
}

timeStamp("Compare each cell subtypes for IPF vs. Control.")
t.disease.list <- sapply(
    simplify=FALSE,
    unique(anno$category),
    function(cat){
        timeStamp(paste("T-test run for", cat))
        # subset by category.
        d <- counts.logcpm[, anno$category == cat]
        a <- anno[anno$category == cat, ]
        # compare IPF vs. Control.
        t <- t_test(d, a$Disease_Identity == "IPF")
        print(head(t))
        return(t)
    }
)



lapply(
    names(t.list),
    function(n){
        d <- t.list[[n]]
        write.csv(d, paste0("tTest_", n, ".csv"))
        return(NULL)
    }
)

save(t.list, file="RData/tTest.RData")

timeStamp("Done!")