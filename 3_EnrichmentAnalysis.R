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

pdf_width <- 11
pdf_height <- 8.5
pdf("go_graph.pdf", width=pdf_width, height=pdf_height)

lapply(contrast_list, function(list){
  lapply(list, function(sublist){
    gost.res <- gost(rownames(sublist), 
                     organism = "hsapiens", 
                     correction_method = "fdr")
    
    gostyplot <- gostplot(gost.res, interactive = F, capped = F)
    
    publish_gostplot(
      gostyplot,
      filename= NULL,
      width = pdf_width,
      height = pdf_height)
    
  })
})
dev.off()

# get gsea for homo sapiens
hs_gsea <- msigdbr(species = "Homo sapiens")

hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)

# C2 - curated gene sets from online pathway databases, publications in PubMed,
#      and knowledge of domain experts.
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", 
                      category = "C2") %>% # msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) 

# C5 -	ontology gene sets consist of genes annotated by the same ontology term
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5") %>% # msigdb collection of interest
  dplyr::select(gs_name, gene_symbol)

# GSEA data input
myDataGSEA <- contrast_list$wt_4h$pos$log2FoldChange %>% 
  'names<-'(rownames(contrast_list$wt_4h$pos)) %>% 
  sort(decreasing = T)

# GSEA analysis
myGSEA.res.c2 <- GSEA(myDataGSEA, TERM2GENE = hs_gsea_c2,verbose = F)
myGSEA.df.c2 <- as_tibble(myGSEA.res.C2@result)

myGSEA.res.c5 <- GSEA(myDataGSEA, TERM2GENE = hs_gsea_c5,verbose = F)
myGSEA.df.c5 <- as_tibble(myGSEA.res.c5@result)

# Visualize GSEA output - interactive table : C2 - curated gene sets
datatable(myGSEA.df.c2, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)
 

gseaplot2(myGSEA.res.c2, 
          geneSetID = 1, 
          pvalue_table = F, 
          title = myGSEA.res.c2$Description[1]) 

# Visualize GSEA output - interactive table : C5 -	ontology gene sets
datatable(myGSEA.df.c5, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

gseaplot2(myGSEA.res.c5, 
          geneSetID = 1, 
          pvalue_table = F, 
          title = myGSEA.res.c5$Description[1]) 


myGSEA.df <- myGSEA.df.c2 %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
