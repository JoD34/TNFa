library(Biobase)
library(clusterProfiler)
library(DT)
library(enrichplot)
library(limma)
library(gplots)
library(gprofiler2)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(tidyverse)

pdf("go_graph.pdf")

lapply(contrasts.list,function(main_list){
  lapply((main_list), function(sub_list){
    gost.res <- gost(rownames(sub_list), 
                     organism = "hsapiens", 
                     correction_method = "fdr")
    gogo <- gostplot(gost.res, interactive = F, capped = F) 
    publish_gostplot(
      gogo, 
      highlight_terms = NA,
      filename = "go_graph.pdf",
      width = NA,
      height = NA)
    print(gogo)
  })
})
dev.off()
