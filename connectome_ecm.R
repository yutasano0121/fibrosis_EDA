library(igraph)
library(RCy3)

sources("functions.R")

timeStamp("Data loading and preprocessing...")

# genes upregulated in IPF
deg <- read.csv("DE/tTest_IPFvsCtrl_combined.csv")
deg <- deg[deg$logFC >= 1, ]  # more than 2-fold upregulation

# ECM genes
ecm <- read.csv("ECM_matrisome_hs_masterlist.csv", stringsAsFactors=FALSE)
ecm <- ecm$Gene.Symbol

# ligand-receptor list
ligRec <- read.csv("PairsLigRec.tsv", sep="\t", stringsAsFactors=FALSE)

# remove rows with a negative pair evidence
ligRec <- ligRec[grep("EXCLUDED", ligRec$Pair.Evidence, invert=TRUE), ]

# remove interactions not related to ECM
lig_ecm <- ligRec$Ligand.ApprovedSymbol %in% ecm
rec_ecm <- ligRec$Receptor.ApprovedSymbol %in% ecm
ligRec <- ligRec[lig_ecm | rec_ecm, ]

# fetch DEGs of Mf, alveolar Mf, and aberrant basaloid
cellType_ofInterest <- c("Macrophage", "Macrophage_Alveolar", "Aberrant_Basaloid")
graph.list <- lapply(
    cellType_ofInterest,
    function(type){
        deg.target <- deg$gene[deg$cellType == type]  # returns a character vector
        deg.others <- deg[deg$cellType != type, ]  # returns a data frame

        # ligands and receptors in Mf-upregulated genes
        deg.target.lig <- deg.target[deg.target %in% ligRec$Ligand.ApprovedSymbol]
        deg.target.rec <- deg.target[scBcells %in% ligRec$Receptor.ApprovedSymbol]

        # Subset the ligand-receptor list with those present in scBcells.
        lig_inTarget <- ligRec$Ligand.ApprovedSymbol %in% deg.target.lig
        rec_inTarget <- ligRec$Receptor.ApprovedSymbol %in% deg.target.rec
        ligRec.sub <- ligRec[lig_inMF | rec_inTarget, ]

        # ligands and receptors present in both Mf and other cells.
        deg.others.lig <- deg.others$gene %in% ligRec.sub$Ligand.ApprovedSymbol
        deg.others.rec <- deg.others$gene %in% ligRec.sub$Receptor.ApprovedSymbol
        deg.others.ligRec <- deg.others[deg.others.lig | deg.others.rec, ]

        target.ligRec.list <- lapply(
            c(deg.target.lig, deg.target.rec),
            function(x){return(c("Macrophage", x))}
        )

        others.ligRec.list <- lapply(
            1:dim(deg.others.ligRec)[1],
            function(i){
                return(c(deg.others.ligRec$cellType[i], deg.others.ligRec$gene[i]))
            }
        )

        # Detect ligand-receptor pairs present between B cells and tissues.
        all.lig <- unique(c(deg.target.lig, deg.others$gene[deg.others.lig]))
        all.rec <- unique(c(deg.target.rec, deg.others$gene[deg.others.rec]))

        fetch.ligRecPairs <- 

        # CAUTION: produces a lot of null values
        ligRecPairs <- lapply(
            1:dim(ligRec.sub)[1], 
            function(i){
                if (
                    (ligRec.sub$Ligand.ApprovedSymbol[i] %in% all.lig) &
                    (ligRec.sub$Receptor.ApprovedSymbol[i] %in% all.rec)
                ){
                    return(c(
                            as.character(ligRec.sub$Ligand.ApprovedSymbol[i]),
                            as.character(ligRec.sub$Receptor.ApprovedSymbol[i])
                    ))
                }
            }
        )


        # make a data frame to make igraph object
        df <- rbind(
            do.call(rbind, target.ligRec.list),
            do.call(rbind, others.ligRec.list),
            do.call(rbind, ligRecPairs)
        )

        # make an igraph object
        g <- graph.data.frame(df)
        V(g)$ligRecPairs <- TRUE
        V(g)$ligRecPairs[1 : length(target.ligRec) + length(others.ligRec)] <- FALSE
        V(g)$category <- "Source"
        V(g)$category[V(g)$name %in% all.lig] <- "Ligand"
        V(g)$category[V(g)$name %in% all.rec] <- "Receptor"
        V(g)$Paired <- FALSE
        V(g)$Paired[V(g)$name %in% as.character(do.call(rbind, ligRecPairs))] <- TRUE

        return(g)
    }
)

# open Cytoscape
cytoscapePing()
createNetworkFromIgraph(g, "cell-tissue connectome")
