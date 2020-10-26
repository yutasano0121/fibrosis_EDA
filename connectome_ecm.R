library(igraph)
library(RCy3)

source("functions.R")

timeStamp("Data loading and preprocessing...")

anno <- read.csv("annotation_QCed_1000gene10MT_reduced.csv", stringsAsFactors=FALSE)
cellTypes <- sapply(
    simplify=FALSE,
    unique(anno$Manuscript_Identity),
    function(type){
        unique(anno$CellType_Category[anno$Manuscript_Identity == type])
    }
)

# genes upregulated in IPF
deg <- read.csv("DE/tTest_IPFvsCtrl_combined.csv", stringsAsFactors=FALSE)
deg <- deg[deg$logFC >= log2(2), ]  # more than 1.5-fold upregulation

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
ligRec_ecm <- ligRec[lig_ecm | rec_ecm, ]

# fetch DEGs related to ligand-receptor ECM interaction
deg_lig <- deg$gene %in% ligRec$Ligand.ApprovedSymbol
deg_rec <- deg$gene %in% ligRec$Receptor.ApprovedSymbol
deg_ligRec <- deg[deg_lig | deg_rec, ]  

# fetch paired ligand and receptor genes present in deg
deg_lig_names <- unique(deg_ligRec$gene[deg_ligRec$gene %in% ligRec$Ligand.ApprovedSymbol])
deg_rec_names <- unique(deg_ligRec$gene[deg_ligRec$gene %in% ligRec$Receptor.ApprovedSymbol])
ligRecPaired <- lapply(  # CAUTION: produces a lot of null values
    1:dim(ligRec_ecm)[1], 
    function(i){
        if (
            (ligRec_ecm$Ligand.ApprovedSymbol[i] %in% deg_lig_names) &
            (ligRec_ecm$Receptor.ApprovedSymbol[i] %in% deg_rec_names)
        ){
            return(c(
                    ligRec_ecm$Ligand.ApprovedSymbol[i],
                    ligRec_ecm$Receptor.ApprovedSymbol[i]
            ))
        }
    }
)
ligRecPaired <- as.data.frame(do.call(rbind, ligRecPaired))

# remove pairs not related to Mf or alveolar Mf
genes_mf <- deg_ligRec$gene[(deg_ligRec$cellType=="Macrophage") | (deg_ligRec$cellType=="Macrophage_Alveolar")]
ligRecPaired <- ligRecPaired[(ligRecPaired$V1 %in% genes_mf) | (ligRecPaired$V2 %in% genes_mf), ]

# remove genes not paired
deg_lig_paired <- deg_ligRec$gene %in% ligRecPaired$V1
deg_rec_paired <- deg_ligRec$gene %in% ligRecPaired$V2
deg_ligRecPaired <- deg_ligRec[deg_lig_paired | deg_rec_paired, ]

# make a category-ligRec pair list for making a graph
deg_ligRecPaired.list <- lapply(
    1:dim(deg_ligRecPaired)[1],
    function(i){
        return(c(deg_ligRecPaired$cellType[i], deg_ligRecPaired$gene[i]))
    }
)
deg_ligRecPaired.list <- as.data.frame(do.call(rbind, deg_ligRecPaired.list))


# make a data frame to make igraph object
df <- rbind(ligRecPaired, deg_ligRecPaired.list)

# make an igraph object
g <- graph.data.frame(df)
V(g)$ligRecPaired <- FALSE  # if each node is a ligand-receptor pair
V(g)$ligRecPaired[1 : length(ligRecPaired)] <- TRUE
V(g)$broadCellType <- "None"  # label each node by broad cell types
for (i in 1:length(V(g))){
    if (V(g)$name[i] %in% names(cellTypes)){
        V(g)$broadCellType[i] <- cellTypes[[V(g)$name[i]]]
    }
}
V(g)$category <- "Source"
V(g)$category[V(g)$name %in% ligRecPaired$V1] <- "Ligand"
V(g)$category[V(g)$name %in% ligRecPaired$V2] <- "Receptor"
#V(g)$Paired <- FALSE
#V(g)$Paired[V(g)$name %in% as.character(do.call(rbind, ligRecPaired))] <- TRUE


# open Cytoscape
cytoscapePing()
createNetworkFromIgraph(g, "ECM connectome")
