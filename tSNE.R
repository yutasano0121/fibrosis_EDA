library(Matrix)
library(edgeR)
library(Rtsne)
library(ggplot2)
library(viridis)

source("function.R")

counts <- readMM(file="counts_QCed_1000gene10MT_reduced.mtx.gz")
anno <- read.csv("annotation_QCed_1000gene10MT_reduced.csv")
ecm <- read.csv("ECM_matrisome_hs_masterlist.csv")

counts <- as.data.frame(log2(cpm(counts) + 1))
counts.ecm <- counts[ecm$Gene.Symbol, ]

set.seed(0)
tsne <- Rtsne(counts.cpm)

p <- plot.tsne(DF=counts.ecm, ANNO=anno)

pngStore(p, "tsne.png")

