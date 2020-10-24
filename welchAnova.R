removeRareGenes <- FALSE
gene_freq_cutoff <- 0.05  # cutoff to decide how many genes are tested

source("functions.R")
timeStamp("RData loading...")
load("RData/ecm_tsne_reduced.RData")
timeStamp("RData loaded.")


# remove PNEC, since it has only 1 observation.
notPNEC <- anno$Manuscript_Identity != "PNEC"
counts.ecm <- counts.ecm[notPNEC, ]
anno <- anno[notPNEC, ]


# divide macrophages by disease states
category <- anno$Manuscript_Identity
Mf_d <- (category == "Macrophage") & (anno$Disease_Identity == "IPF")
Mf_n <- (category == "Macrophage") & (anno$Disease_Identity == "Control")
MfA_d <- (category == "Macrophage_Alveolar") & (anno$Disease_Identity == "IPF")
MfA_n <- (category == "Macrophage_Alveolar") & (anno$Disease_Identity == "Control")

category[Mf_d] <- "Macrophage_IPF"
category[Mf_n] <- "Macrophage_Control"
category[MfA_d] <- "Mf_Alveolar_IPF"
category[MfA_n] <- "Mf_Alveolar_Control"

anno$category <- category


test_anova <- function(
    counts=counts.ecm, meta=anno,
    saveResult=TRUE, output_name=""
){
    if (removeRareGenes){
        # remove genes expressed by less than 5 % of all cells
        expressed_genes <- apply(
            counts,
            1,
            function(x){
                sum(x != 0) / dim(counts)[2] >= gene_freq_cutoff
            }
        )
        counts <- counts[expressed_genes, ]
    }

    test.result <- sapply(
        colnames(counts),
        function(gene){
            d <- data.frame(as.numeric(counts[, gene]), meta$category)
            colnames(d) <- c("gene", "group")
            test <- oneway.test(gene~group, data=d)
            return(test$p.value)
        }
    )
    test.result.sub <- test.result[!is.na(test.result)]
    fdr <- p.adjust(test.result.sub, "BH")

    df <- data.frame(fdr, test.result.sub)
    colnames(df) <- c("fdr", "p")
    df <- df[order(df$fdr), ]
    print(
        paste(
            sum(df$fdr <= 0.05), "DEGs found out of", dim(df)[1], "genes tested!"
        )
    )

    write.csv(df, output_name)
    return(df)
}


# t-test to identify cluster markers
t_test <- function(DF, CATEGORY){
    t.result <- sapply(
        colnames(DF),
        function(gene){
            t <- t.test(DF[, gene]~CATEGORY)
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
    t.result <- t.result[(t.result$fdr <= 0.05) & (t.result$logFC > 0), ]

    return(t.result)
}

t.list <- sapply(
    simplify=FALSE,
    unique(anno$category),
    function(cat){
        timeStamp(paste("T-test run for", cat))
        t <- t_test(counts.ecm, anno$category==cat)
        print(head(t))
        return(t)
    }
)
timeStamp("T-test run for Macrophages")
t.list$Mf <- t_test(counts.ecm, anno$Manuscript_Identity == "Macrophage")


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