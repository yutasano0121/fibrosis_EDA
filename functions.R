library(ggplot2)
library(viridis)

timeStamp <- function(comment)
{
    ts <- paste(comment, Sys.time())
    print(ts)
}


pngStore <- function(plot_data, output_name, HEIGHT, WIDTH)
{
    png(output_name, height = HEIGHT, width = WIDTH, unit = "in", res = 600)
    print(plot_data)
    dev.off()
}


# find genes expressed at a frequency higher than the threshold
find_expressedGenes <- function(DATA, FREQUENCY)  # data needs to be genes x cells
{
    gene.count <- rowSums(DATA != 0)
    detection.rate <- gene.count / dim(DATA)[2]
    expressed_genes <- detection.rate > FREQUENCY  # default 0.1
    return(expressed_genes)
}

# select high-variance genes in the expression matrix (genes x cells)
subset_genes <- function(data, geneNum)
{
    gene_var <- apply(data, 1, var)
    data_ordered <- data[order(gene_var, decreasing=TRUE), ]
    data_ordered <- data_ordered[1:geneNum, ]
    return(data_ordered)
}


# find most highly expressed Ig isotype. pass it to apply(counts, 2, findIg)
findIg <- function(x){
    ensg <- HC[which.max(x)]
    return(anno.geneNames$gene_symbol[anno.geneNames$ENSG==ensg])
}


# a function to extract differentially expressed gene names
extract_DEGs <- function(TAG, TOP=NULL){
    if (is.null(TOP)){
        deg <- rownames(TAG)[TAG$threshold]
    } else {  # take only {TOP} genes with lowest FDR
        tag <- TAG[TAG$threshold, ]
        tag <- tag[order(tag$FDR), ]
        deg1 <- rownames(tag)[tag$logFC > 0][1:TOP]
        deg2 <- rownames(tag)[tag$logFC < 0][1:TOP]
        deg <- c(deg1, deg2)
    }

    return(deg)
}


# Hierarchical-cluster DEGs using a selected method.
# Returns
hierarchicalCluster <- function(counts, clusteringMethod, n_cluster){
    # Counts = row: genes, col: groups or cells
    hc_gene <- hclust(as.dist(1 - cor(t(counts))),  method = clusteringMethod)
    # Divid the tree into clusters.
    mycl <- cutree(hc_gene, k=n_cluster)
    # make it a data frame to feed to pheatmap function
    mycl <- as.data.frame(as.numeric(mycl))
    rownames(mycl) <- rownames(counts)
    # make another data frame to be saved as .RData
    geneClusters <- mycl
    geneClusters$symbol <- rownames(mycl)
    colnames(mycl) <- c("cluster")
    colnames(geneClusters) <- c("cluster", "symbol")
    geneClusters <- geneClusters[hc_gene$order, ]  # reorder it to the order of heatmap
    # save geneClusters and the hierachical clustering info as .RData

    return(
        list(
            geneClusters=geneClusters,
            mycl=mycl,
            hc_gene=hc_gene
        )
    )
}


# Plot a heatmap of DEG clusters.
# Color break needs to be defined by 'quantile_breaks' function.
plot_heatmap <- function(
    counts, hc_gene, mycl, n_cluster, col_break,
    rowClusterAnnotation=TRUE  # if cluster annotation is added to the heatmap.
){
    if (rowClusterAnnotation){
        pheatmap(
            counts,
            color=plasma(length(col_break) - 1),
            breaks=col_break,
            cluster_rows=hc_gene,
            cluster_cols=FALSE,
            show_rownames=FALSE,
            angle_col=45,
            fontsize=14,
            annotation_row=mycl,
            annotation_colors=list(cluster=grey.colors(n_cluster)),
            annotation_legend=FALSE,
            annotation_names_row=FALSE
        )
    }
    else {
        pheatmap(
            counts,
            color=plasma(length(col_break) - 1),
            breaks=col_break,
            cluster_rows=hc_gene,
            cluster_cols=FALSE,
            show_rownames=FALSE,
            angle_col=45,
            fontsize=14,
            annotation_legend=FALSE,
            annotation_names_row=FALSE
        )
    }
}


makeRownamesENTREZ <- function(DF, ORG.DB="org.Hs.eg.db")
{
    ids <- bitr(rownames(DF), fromType="SYMBOL", toType="ENTREZID", OrgDb=ORG.DB)
    # fetch duplicated gene symbols
    if (sum(duplicated(ids$SYMBOL)) > 0 ){
        duplicated.symbols <- ids$SYMBOL[duplicated(ids$SYMBOL)]
        print(paste("Gene symbols with multiple ENTREZID are removed:", paste(duplicated.symbols)))
        # remove them
        ids <- ids[!ids$SYMBOL %in% duplicated.symbols, ]
    }
    DF.out <- DF[rownames(DF) %in% ids$SYMBOL, ]
    # change rownames to ENTREZID
    ids <- ids[order(ids$SYMBOL), ]
    DF.out <- DF.out[order(rownames(DF.out)), ]
    rownames(DF.out) <- ids$ENTREZID

    return(DF.out)
}


# run GO-enrichment test on the list of differentially expressed genes
test.GOenrichment <- function(GENES,
                              ONT='BP', P.ADJUST='BH',
                              ORG.DB='org.Hs.eg.db',  # 'Hs' or 'Mm'
                              P.CUTOFF=0.05, Q.CUTOFF=0.2,
                              MIN.SIZE=10, MAX.SIZE=500,
                              READABLE=TRUE)
{
    ego <- enrichGO(gene          = GENES,
                    OrgDb         = ORG.DB,
                    keyType       = 'ENTREZID',
                    ont           = ONT,
                    pAdjustMethod = P.ADJUST,
                    pvalueCutoff = P.CUTOFF,
                    qvalueCutoff  = Q.CUTOFF,
                    minGSSize     = MIN.SIZE,
                    maxGSSize     = MAX.SIZE,
                    readable      = READABLE)
    return(ego)
}

# run GO-enrichment test on the list of differentially expressed genes
test.KEGGenrichment <- function(GENES, P.ADJUST='BH',
                                ORG.DB='hsa',  # 'hsa' or 'mmu'
                                MIN.SIZE=10, MAX.SIZE=500,
                                P.CUTOFF=0.05)
{
    ekegg <- enrichKEGG(gene          = GENES,
                      organism      = ORG.DB,
                      pAdjustMethod = P.ADJUST,
                      pvalueCutoff  = P.CUTOFF,
                      minGSSize     = MIN.SIZE,
                      maxGSSize     = MAX.SIZE)
    return(ekegg)
}



# run go-GSEA (input needs to be a ranked gene lists)
test.GOgse <- function(rankedGENES,
                      ONT='BP', P.ADJUST='BH',
                      P.CUTOFF=0.05, PERM = 1000,
                      MIN.SIZE=100, MAX.SIZE=500)
{
    gsego <- gseGO(geneList    = rankedGENES,
                  OrgDb        = 'org.Hs.eg.db',
                  keyType      = 'ENTREZID',
                  ont          = ONT,
                  pAdjustMethod = P.ADJUST,
                  nPerm        = PERM,
                  minGSSize    = MIN.SIZE,
                  maxGSSize    = MAX.SIZE,
                  pvalueCutoff = P.CUTOFF)
    return(gsego)
}


# Fetch genes of a GO/KEGG pathway of interest.
fetchGenes <- function(
    GROUP,  # groupGO  or enrichKEGG result
    PATHWAY,  # pathway name
    KEGG=FALSE  # if input is an enrichKEGG result
){
    # fetch a result data frame from a GO/KEGG output
    df <- as.data.frame(GROUP@result)
    # extract genes assigned to 'innate immune response' GO term
    genes <- df$geneID[df$Description == PATHWAY]
    genes <- unlist(strsplit(genes, "/"))
    if (KEGG){
        ids <- bitr(
            genes,
            fromType="ENTREZID",
            toType="SYMBOL",
            OrgDb="org.Hs.eg.db"
        )
        genes <- ids$SYMBOL
    }
    return(genes)
}


# run KEGG-GSEA (input needs to be a ranked gene lists)
test.KEGGgse <- function(rankedGENES,
                      ONT='BP', P.ADJUST='BH',
                      P.CUTOFF=0.05, PERM = 1000,
                      MIN.SIZE=100, MAX.SIZE=500)
{
    gsekegg <- gseKEGG(gene          = rankedGENES,
                      organism      = 'hsa',
                      pAdjustMethod = P.ADJUST,
                      pvalueCutoff  = P.CUTOFF,
                      nPerm         = PERM,
                      minGSSize     = MIN.SIZE,
                      maxGSSize     = MAX.SIZE)
    return(gsekegg)
}


# calculate signal-to-noise ratio
stn <- function(VECTOR, CATEGORY)
{
    comparison <- levels(CATEGORY)
    sd.1 <- sd(VECTOR[CATEGORY==comparison[1]])
    sd.2 <- sd(VECTOR[CATEGORY==comparison[2]])
    mean.1 <- mean(VECTOR[CATEGORY==comparison[1]])
    mean.2 <- mean(VECTOR[CATEGORY==comparison[2]])

    stnRatio <- (mean.2 - mean.1) / (sd.1 + sd.2)
    return(stnRatio)
}

# calculate Baumgartner-Weiss-Schindler statistics
bws <- function(VECTOR, CATEGORY, METHOD="BWS")  # METHOD = 'BWS' or 'NEuhauser'
{
    comparison <- levels(CATEGORY)
    vec1 <- VECTOR[CATEGORY==comparison[1]]
    vec2 <- VECTOR[CATEGORY==comparison[2]]
    # one-sided BWS test
    out <- bws_test(vec1, vec2, method=METHOD)
    if ((METHOD == 'BWS') & (mean(vec2) - mean(vec1) < 0)){out$statistic <- -out$statistic}
    return(out$statistic)
}


# set color breaks for heatmap according to quantile distribution
quantile_breaks <- function(DATA, n = 100){
    breaks <- quantile(DATA, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
}


# run an over-representation (hypergeometric) test
test.enrichment <- function(GENES, TARGET.LIST=likeAhnak, BACKGROUND.SIZE=geneNum.immgen, TARGET.SIZE=NULL)
{
    geneNum.deg <- length(GENES)
    numHits <- sum(GENES %in% TARGET.LIST)
    if (is.null(TARGET.SIZE)){
        num.target <- length(TARGET.LIST)
    } else {num.target <- TARGET.SIZE}
    p.value <-  phyper(q=numHits, k=geneNum.deg,
                       m=num.target, n=BACKGROUND.SIZE-num.target,
                       lower.tail=FALSE)
    return(list(p.value=p.value, ratio=numHits / geneNum.deg, size=geneNum.deg))
}


# fetch the default color palette of ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# volcano plot
# take an edgeR result as an input
plot.volcano <- function(TAG, num_label_genes=10)
{
    df <- as.data.frame(TAG)
    p <- ggplot(data = df, aes(x = logFC, y = -log10(FDR), colour = threshold))
    p <- p + geom_point(alpha = 0.4, size = 1.75)
    p <- p + theme(legend.position = "none")
    p <- p + xlab(paste0("log2 fold change (Kidney / Tonsil)"))
    p <- p + ylab("-log10 FDR")
    #p <- p + xlim(c(-10, 10))
    p <- p + scale_colour_manual(values = c("darkgrey", "red"))

    # label top hits (arbitrary)
    #tmp <- df[order(-abs(df$logFC)), ]  # order genes by logFC
    tmp <- df[order(df$FDR), ]  # order genes by fdr
    #tmp1 <- rownames(tmp[tmp$logFC > 0, ][1:4, ])
    #tmp2 <- rownames(tmp[tmp$logFC < -7.5, ])
    tmp1 <- rownames(tmp[tmp$logFC > 0, ][1:num_label_genes, ])
    tmp2 <- rownames(tmp[tmp$logFC < 0, ][1:num_label_genes, ])
    #tmp3 <- rownames(tmp[tmp$logFC < 0, ][1:10, ])
    #tmp4 <- rownames(tmp[tmp$logFC > 0, ][1:7, ])
    #tmp <- c(tmp1, tmp2, tmp3, tmp4)
    tmp <- c(tmp1, tmp2)
    tmp <- tmp[!duplicated(tmp)]
    p <- p + ggrepel::geom_text_repel(
        data=df[tmp ,],
        aes(x=logFC, y=-log10(FDR), label=tmp),
        segment.color="darkgray",
        colour='black',
        force=10
    )  # don't forget x and y aesthetics variables here!

    # further format labeling
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.position="none")


    return(p)
}


# simply label top 10 hits (both FDR and logFC, so max 20 labels per side)
plot.volcano2 <- function(TAG)
{
    df <- as.data.frame(TAG)
    p <- ggplot(data = df, aes(x = logFC, y = -log10(FDR), colour = threshold))
    p <- p + geom_point(alpha = 0.4, size = 1.75)
    p <- p + theme(legend.position = "none")
    p <- p + xlab(paste0("log2 fold change (Kidney / Tonsil)"))
    p <- p + ylab("-log10 FDR")
    #p <- p + xlim(c(-10, 10))
    p <- p + scale_colour_manual(values = c("darkgrey", "red"))

    # label top hits (top10 FDR and logFC)
    #tmp <- df[order(-abs(df$logFC)), ]  # order genes by logFC
    tmp <- df[order(df$FDR), ]  # order genes by fdr
    tmp1 <- rownames(tmp[tmp$logFC > 0, ][1:10, ])
    tmp2 <- rownames(tmp[tmp$logFC < 0, ][1:10, ])
    tmp <- df[order(-abs(df$logFC)), ]
    tmp3 <- rownames(tmp[tmp$logFC > 0, ][1:10, ])
    tmp4 <- rownames(tmp[tmp$logFC < 0, ][1:10, ])
    tmp <- c(tmp1, tmp2, tmp3, tmp4)
    tmp <- tmp[!duplicated(tmp)]
    p <- p + ggrepel::geom_text_repel(data=tag[tmp ,], aes(x=logFC, y=-log10(FDR), label=tmp), segment.color = "darkgray", colour='black', force=10)  # don't forget x and y aesthetics variables here!

    # further format labeling
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.position="none")


    return(p)
}

# plot enriched GO/KEGG terms
# take an enrichGO/enrichKEGG result as an input
plot.enrichedTerms <- function(eRESULT, showCategory=NULL)
{
    df <- eRESULT@result
    df <- df[df$p.adjust < 0.05, ]  # subset significant terms
    df <- df[order(df$p.adjust), ]  # order the data frame by FDR
    if (!is.null(showCategory)){
        num_significant_terms = dim(df)[1]
        # Subset rows by a smaller value of showCategory or number of enriched terms.
        df <- df[1:min(num_significant_terms, showCategory), ]}
        # Round p.adjust values.
        df$p.adjust <- signif(df$p.adjust, digit = 2)
    p <- ggplot(data = df, aes(y = Count, x = reorder(Description, -p.adjust), fill = p.adjust))
    p <- p + geom_bar(stat = "identity") + coord_flip()
    p <- p + scale_fill_viridis(direction=-1, option="plasma")

    # format labeling
    p <- p + xlab("") + ylab("Gene count") + labs(fill = "FDR")
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14))

    return(p)
}


# themes for scatter plots
theme.scatter <- theme_classic() +
                 theme(axis.title = element_text(size=14, colour='black'),
                       axis.text = element_text(size=12, colour="black"),
                       legend.title = element_text(size = 14),
                       legend.text = element_text(size = 14))


# make a t-SNE plot
# change .LAB parameters as well when aesthetics parameters are changed
plot.tsne <- function(
    DF, ANNO,
    X='V1', Y='V2',
    COLOR='CellType_Category', COLOR.LAB='Cell types',
    GENE=FALSE, GENE.DF=None
)
{
    p <- ggplot(DF, aes(x=DF[, X], y=DF[, Y]))
    p <- p + theme.scatter
    p <- p + xlab("tSNE-1") + ylab("tSNE-2")

    if (GENE){
        p <- p + geom_point(aes(color=GENE.DF[, COLOR]))
        p <- p + scale_color_viridis_c()
        if (COLOR.LAB=="Patient"){
            # default
            p <- p + labs(colour=COLOR)
            p <- p + theme(legend.title = element_text(face="italic"))
        } else {p <- p + labs(colour=COLOR.LAB)}
    } else{
        p <- p + geom_point(aes(color=ANNO[, COLOR]))
        p <- p + labs(colour=COLOR.LAB)
    }

    return(p)
}


plot.violin <- function(
    DF, ANNO, Y, X='sample',
    FILL=NULL, FILL.LAB="Add a label", JITTER=FALSE, GENE=TRUE, Y.LAB=NULL)
{
    if (is.null(FILL)){
        p <- ggplot(DF, aes(x=ANNO[, X], y=DF[, Y]))
    } else {
        p <- ggplot(DF, aes(x=ANNO[, X], y=DF[, Y], fill=ANNO[, FILL]))
        p <- p + labs(fill = FILL.LAB)
    }

    p <- p + geom_violin(scale="width")

    if (JITTER){
        if (is.null(FILL)){
            p <- p + geom_point(size=0.5, alpha=0.7, position=position_jitter(width=0.35))
        } else{
            p <- p + geom_point(size=0.5, alpha=0.7, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.7))
        }
    } else {p <- p + geom_boxplot(width=0.1, position=position_dodge(0.9))}


    # format labeling
    p <- p + theme_classic()
    p <- p + theme(
        axis.title.x = element_text(size=14, colour='black'),
        axis.title.y = element_text(size=14, colour='black'),
        axis.text.x = element_text(size = 12, colour="black", angle=45, hjust=1),
        axis.text.y = element_text(size = 12, colour="black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )
    p <- p + xlab("") + ylab(Y)

    # make gene names italic
    if (GENE==TRUE){
        p <- p + theme(axis.title.y = element_text(size=14, colour='black', face="italic"))
    } else {p <- p + ylab(Y.LAB)}

    return(p)
}

# make a bar plot with error bars
plot.bar <- function(DF, X, Y='mean', X.IN.DF=TRUE, ERRORBAR=TRUE, ERROR='se', Y.LAB="Score", WIDTH=NULL)
{
    if (X.IN.DF == FALSE)
    {
        # default. X-axis is determined by another data frame, list, or vector
        p <- ggplot(DF, aes(x=X, y=DF[, Y])) + geom_bar(stat="identity", width=WIDTH)
    } else {
        # if DF contains a column of a categorical variable to determine the x-axis
        p <- ggplot(DF, aes(x=DF[, X], y=DF[, Y])) + geom_bar(stat="identity", width=WIDTH)
    }

    # add error bars if needed
    if (ERRORBAR == TRUE)
    {
        p <- p + geom_errorbar(aes(ymin=DF[, Y] - DF[, ERROR], ymax=DF[, Y] + DF[, ERROR]),
                                    width=.2, position=position_dodge(.9))
    }

    # format labeling
    p <- p + xlab("") + ylab(Y.LAB)
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size=14, colour='black'),
                axis.text.x = element_text(size = 14, colour="black", angle=45, hjust=1),
                axis.text.y = element_text(size = 12, colour="black"))

    return(p)
}


plot.enrichNetwork <- function(enrichResult,
                                showCategory = 10,
                                kegg=FALSE,
                                saveCytoscapeObj=FALSE,
                                nameCytoscapeObj="temp")
{
    # make a data frame from an enrichResult object
    df <- as.data.frame(enrichResult)[1:showCategory, ]

    # make a list of genes for top enriched terms
    termsAndGenes <- strsplit(df$geneID, "/")
    names(termsAndGenes) <-  df$Description

    # make the list into a data frame for making a graph
    df_termsAndGenes <- lapply(names(termsAndGenes), function(name){
        data.frame(categoryID=rep(name, length(termsAndGenes[[name]])),
                   Gene=termsAndGenes[[name]])
        })
    df_termsAndGenes <- do.call('rbind', df_termsAndGenes)

    if (kegg == TRUE)
    {
        # change geneID from ENTREZID to SYMBOL
        library(org.Hs.eg.db)

        ids <- bitr(geneID=df_termsAndGenes$Gene,
                    fromType="ENTREZID", toType="SYMBOL",
                    OrgDb="org.Hs.eg.db",
                    drop=FALSE)
        df_termsAndGenes$Gene <- sapply(df_termsAndGenes$Gene, function(x){return(ids$SYMBOL[ids$ENTREZID==x])})
    }

    # convert the data frame into a graph
    g <- graph.data.frame(df_termsAndGenes, directed=FALSE)

    V(g)$size <- 1
    V(g)$size[1:showCategory] <- 4
    fdr <- df$p.adjust
    V(g)$color <- 0
    V(g)$color[1:showCategory] <- as.numeric(fdr)
    #my.palette <- viridis(n=100, begin=0.2, end=0.9, direction=-1, option="plasma")
    my.palette <- viridis(n=100, direction=-1, option="plasma")


    if(saveCytoscapeObj==TRUE)
    {
        gg <- g
        V(gg)$category <- FALSE
        V(gg)$category[1:showCategory] <- TRUE
        #fdr.scaled <- ((fdr - min(fdr)) / (max(fdr) - min(fdr))) * 99 + 1  # scale between 1 and 100
        if (kegg){
            #V(gg)$color[showCategory + 1] <- 0.0429
           # print ("color scale max is set as 0.0429")
        }
        library(RCy3)
        cytoscapePing()
        createNetworkFromIgraph(gg, nameCytoscapeObj)
        timeStamp("igraph transfered to Cytoscape.")
    }

    V(g)$name[1:showCategory] <-  paste0("bold('", df$Description, "')")

    p <- ggraph(g, layout="kk") +
        geom_edge_link(alpha=.5, colour='darkgrey') +
    #    geom_node_point(aes_(color=~color), size=V(g)$size, shape=21) +
        geom_node_point(aes_(color=~color), size=V(g)$size) +
        scale_color_gradientn(name = "FDR", colors=my.palette, na.value = "dimgrey")
    coord <- layout_with_kk(g)
    p <- p + geom_text_repel(x=coord[, 1], y=coord[, 2], label=V(g)$name, segment.color="dimgrey", segment.linetype=3, parse=TRUE, box.padding=0.5, point.padding=0.4, max.iter=10000)
    p <- p + theme_void()

    return(p)
}


scaleColorLimit <- function(LIST.E.OBJECT)
{
    fetch.Plimits <- function(E.OBJECT)
    {
        df <- as.data.frame(E.OBJECT)
        df <- df[df$p.adjust < 0.05, ]  # significant ones
        maxP <- max(df$p.adjust)
        minP <- min(df$p.adjust)  # highest p-value
        return(c(max=maxP, min=minP))
    }

    Plimits <- lapply(LIST.E.OBJECT, fetch.Plimits)
    Plimits <- as.data.frame(do.call(rbind, Plimits))
    upperLimit <- max(Plimits$max)
    lowerLimit <- min(Plimits$min)
    limitList <- c(low=lowerLimit * (1 - 1E-6), high=upperLimit * (1 + 1E-6))

    if (FALSE)
    {
        # a script to make a color scale
        myPalette <- viridis(100, direction=-1, option="plasma")
        scale_fill_gradientn(colors=myPalette, limits=limitList, breaks=c())
    }

    return(limitList)
}


# plot a pie chart from annotation data frame
plot.pie <- function(ANNO, VAR='Ig.class', LEGEND.TITLE="Ig class")
{
    if (class(ANNO[, VAR]) != "factor") {stop("VAR needs to be a factor.")}

    # fetch a summary of the variable of interest
    summary.VAR <- summary(ANNO[, VAR])
    # make a data frame containing count and category info
    y = data.frame(count=summary.VAR, category=factor(names(summary.VAR), levels=levels(ANNO[, VAR])))
    y$count.rev <- rev(y$count)  # for text labeling
    # make a polar bar plot
    bp <- ggplot(y, aes(x=factor(1), y=count, fill=category))
    bp <- bp + geom_bar(stat="identity")
    pie <- bp + coord_polar("y", direction=-1)
    blank_theme <- theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            legend.title=element_text(size=14),
            legend.text=element_text(size=14))
    library(scales)
    p <- pie + blank_theme
    #p <- p + ggrepel::geom_text_repel(aes(y = count.rev/2 + c(0, cumsum(count.rev)[-length(count.rev)]), label = percent(count.rev/sum(count.rev))), size=6, box.padding=1)
    p <- p + guides(fill=guide_legend(title=LEGEND.TITLE))
    #p <- p + scale_fill_manual(values=colorPal)

    return(p)
}

# calculate a scaled sum of expression levels
calcScore <- function(DF, GENES)  # data = genes x cells, genes = character vector. NOT A FACTOR VECTOR!
{
    scaled <- apply(DF[GENES, ], 1, scale)  # scale expression (row-wise)
    # scale() transposes the matrix!
    score <- apply(scaled, 1, sum)  # sum up (row(= cell)-wise))
    return(score)
}
