saveImage <- TRUE
trimSize <- 0  # percentile to trim when calculating mean
clustering <- FALES
num_cluster <- 6

# packages
library(edgeR)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(pheatmap)


# Load functions.
source("functions.R")

# Load lists of marker genes.
deg <- read.csv("DE/tTest_markers_combined.csv", stringsAsFactors=FALSE)

# Load count data (10%)
load("RData/")
print("Data loaded.")


# reorder the annotation
ord <- with(anno, order(Manuscript_Identity, cellType_Category))
anno <- anno[ord, ]

# Take mean values in each category
counts_mean <- aggregate(
    t(counts),
    by=list(anno$Manuscript_Identity),
    FUN=function(x){mean(x, trim=trimSize)}
)
rm(counts)

rownames(counts_mean) <- counts_mean$Group.1
counts_mean <- t(counts_mean[, -1])  # make it gene x category


# Scale the data across cells to have unit variance
counts_scaled <- apply(counts_mean, 1, scale)  # The result is categories x genes
rownames(counts_scaled) <- colnames(counts_mean)
counts_scaled <- t(counts_scaled)  # genes x categories

# set quantile color breaks
col_break <- quantile_breaks(counts_scaled)

print("Data processed.")


print("Make a heatmap.")

if (clustering){
    # make a hierachical cluster of rows
    method_list <- c("complete", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
    #method_list <- c("ward.D2")


    hc_list <- sapply(
        method_list,
        function(method){
            hierarchicalCluster(counts_scaled, method, num_cluster)
        },
        simplify=FALSE
    )

    heatmap_list <- sapply(
        simplify=FALSE,
        method_list,
        function(method){
            if (saveImage){
                png(
                    paste0(imageDir, "/heatmap4comps_", method, "_", num_cluster, "_clusters.png"),
                    height=8, width=4.1, units="in", res=300
                )
                plot_heatmap(
                    counts_scaled, hc_list[[method]]$hc_gene, hc_list[[method]]$mycl,
                    num_cluster, col_break
                )
                dev.off()
            } else {
                plot_heatmap(
                    counts_scaled, hc_list[[method]]$hc_gene, hc_list[[method]]$mycl,
                    num_cluster, col_break
                )
            }
        }
    )

    # Save a list of clustering results.
    save(
        hc_list, num_cluster,
        file=paste0(RDataDir, "/hc_4comps_", num_cluster, "clusters.RData")
    )
} else {
    png(
        "image/heatmap_markers_mean_reduced.png",
        height=8, width=8, units="in", res=300
    )
    pheatmap(
        counts_mean,
        color=viridis(length(col_break) - 1),
        breaks=col_break,
        cluster_rows=FALSE,
        cluster_cols=FALSE,
        show_rownames=FALSE,
        angle_col=45,
        fontsize=14,
    )
    dev.off()
}

print("Done.")
