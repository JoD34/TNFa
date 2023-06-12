# Import libraries ----
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(plotly)
library(pheatmap)
library(RColorBrewer)
library(ReportingTools)
library(scales)
library(stringr)
library(tidyverse)
library(vsn)

# Functions
getNames <- function(list){
  name <- names(list) %>% 
    str_extract( "([0-9]+h|wt)_([0-9]+h|wt)") %>% 
    strsplit("_") %>% 
    unlist() %>% 
    sapply(function(elem){
      elem <- if(!(elem %in% "wt")) 
        paste0("TNF","\u03B1","-",elem) else "Wild Type"
    })
}

# Open pdf to save graphs
pdf("graph_tnfa.pdf", width = 11, height = 8.5)

# Make the design table ----
design.table <- data.frame(
  condition = str_extract(col.names, "TNFa_[0-9]+h|WT")[1:6],
  sample = str_extract(col.names, "(1|2)$")[1:6]
) %>% 
  mutate(
    condition = factor(condition, levels=c("WT", "TNFa_4h", "TNFa_24h")),
    sample = factor(sample)) %>% 
  'rownames<-'(col.names)

# DESeq2 data generation ----
dds <- DESeqDataSetFromTximport(txi.kallisto, 
                                colData = design.table, 
                                design = ~ condition + sample)

# Run DESeq2
dds <- DESeq(dds)

# Implement contrast between conditions + correct for NAs
graph_list <- list(
con_24h_4h <- results(dds,contrast=c("condition", "TNFa_4h", "TNFa_24h")) %>% 
  na.omit(),
con_wt_4h <- results(dds, contrast=c("condition", "TNFa_24h", "WT")) %>% 
  na.omit(),
con_wt_24h <- results(dds, contrast=c("condition", "TNFa_4h", "WT")) %>% 
  na.omit()
)

# Keep  genes with [padj < 0.01] & [log2FC > 1.5]
visualize_data <- list()
for(i in seq_along(graph_list){
  filtered <- graph_list[[i]][con_24h_4h$padj < 0.01 &
                              abs(con_24h_4h$log2FoldChange) > log2(1.5), ]
}

# Order by padj
con_24h_4h.order <- con_24h_4h.filt[order(con_24h_4h.filt$padj), ]
con_wt_4h.order <- con_wt_4h.filt[order(con_wt_4h.filt$padj), ]
con_wt_24h.order <- con_wt_24h.filt[order(con_wt_24h.filt$padj), ]

# Seperate the upregulated from the downregulated
list_24h_4h <- list(
  pos = con_24h_4h.order[con_24h_4h.order$log2FoldChange > 0, ],
  neg = con_24h_4h.order[!con_24h_4h.order$log2FoldChange > 0, ]
)
list_wt_4h <- list(
  pos <- con_wt_4h.order[con_wt_4h.order$log2FoldChange > 0, ],
  neg <- con_wt_4h.order[!con_wt_4h.order$log2FoldChange > 0, ]
)
list_wt_24h <- list(
  pos <- con_wt_24h.order[con_wt_24h.order$log2FoldChange > 0, ],
  neg <- con_wt_24h.order[!con_wt_24h.order$log2FoldChange > 0, ]
)

# Regroup in a list for further usage
contrasts.list <- list(con_24h_4h=list_24h_4h,
                       con_wt_4h=list_wt_4h,
                       con_wt_24h=list_wt_24h)

# Volcano plots ----

for(i in seq_along(graph_list)){
  names <- getNames(graph_list[i])

  value <- graph_list[[i]]
  labels <- data.frame(value[order(log10(value$padj)), ][1:25,])
  
  myPlot <- ggplot(as_tibble(value)) +
    aes(x = log2FoldChange, 
        y = -log10(padj)) +
    geom_point(aes(color = log2FoldChange > 0),
               show.legend = F, 
               shape = 1, 
               size = 2) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = 0.5) +
    labs(title = paste0("Volcano plot: ", name[[2]], " in relation to ", name[[1]]),
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("adjusted pValue")),
         caption = paste("Produced on", Sys.time())) +
    theme_bw()
  
  y <- diff(layer_scales(myPlot)$y$range$range)/2
  
  myPlot <- myPlot +
    annotate("text",
             label = "Downregulated", 
             alpha = 0.5, 
             angle = 90, 
             x = -0.5, 
             y = y,
             size = 3) +
    annotate("text",
             label = "Upregulated", 
             alpha = 0.5, 
             angle = -90, 
             x = 0.5, 
             y = y,
             size = 3)
    print(myPlot)
}
# Apply log fold change shrinkage -> lower noises from log2(FC)----
ma_data <- list(
con_wt_4h.noNA.filt = lfcShrink(dds,
                                coef=c("condition_TNFa_4h_vs_WT"), 
                                type = "apeglm") %>% 
  na.omit(),
con_wt_24h.noNA.filt = lfcShrink(dds,
                                 coef=c("condition_TNFa_24h_vs_WT"), 
                                 type = "apeglm") %>% 
  na.omit()
)
for (i in seq_along(ma_data)){
  names <- getNames(ma_data[i])
  value <- ma_data[[i]]
  ggplot(value)
  myMA <- plotMA(value, colSig= "royalblue4", alpha = 0.01, 
         main=paste0("MAplot: ", names[2], " in relation to ", names[1]),
         ylab="log (FC)") 
    break
}

# Effects of transformations on the variance ----

# vst -> variance stabilizing transformation
vsd <- vst(dds, blind=F)
meanSD <- meanSdPlot(assay(vsd))
print(meanSD)

# PCA analysis
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
namesPCA <- attr(pcaData,"names")

# PCA plot
pca <- ggplot(pcaData) +
  aes(PC1, PC2, color=condition, shape=sample) +
  geom_point(size=3) +
  xlab(paste0(namesPCA[1], ": ", percentVar[1], " %")) +
  ylab(paste0(namesPCA[2], ": ", percentVar[2], " %")) + 
  coord_fixed() +
  theme_bw()
print(pca)

# Heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds))

heat <- pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
print(heat)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Generate heatmap of the distance matrix 
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(name="BuPu", n=9)))(255)
dist_mat <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample-to-sample Distance visualization")
print(dist_mat)

grDevices::dev.off()
