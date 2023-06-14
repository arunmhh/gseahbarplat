library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(tidyverse)

# SET THE DESIRED ORGANISM HERE
#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)
#library(org.Hs.eg.db)

# reading in data from excle.csv
df_sp = read.csv("sp_genelist.csv", header = TRUE)

# we want the log2 fold change 
#sp_gene_list <- df_sp$Log2FoldChange

#name the vector
#names(sp_gene_list) <- df_sp$Gene.symbol

# omit any NA values 
#sp_gene_list<-na.omit(sp_gene_list)

# sort the list in decreasing order 
#(required for clusterProfiler)
#sp_gene_list = sort(sp_gene_list, decreasing = TRUE)
##### gene symblo to Entrez id convesion
gene.df <- bitr(df_sp$Gene.symbol,
     fromType = "SYMBOL",
     toType = c("ENTREZID"),
     OrgDb = org.Hs.eg.db)

gene.df <- rename(gene.df,Gene.symbol= SYMBOL)

## now join the data
df_sp_entrez <- inner_join(df_sp,gene.df, by = "Gene.symbol")

genelist <- df_sp_entrez$Log2FoldChange
genelist1 <- df_sp_entrez$Gene.symbol

names(genelist) <- as.character(df_sp_entrez$ENTREZID)

genelist = sort(genelist, decreasing = T)


####Gene Set Enrichment
gse <- gseGO(geneList =genelist, 
             OrgDb = org.Hs.eg.db,
             ont ="CC", 
             keyType = "ENTREZID", 
             #nPerm = 10000, 
             minGSSize = 100, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = FALSE,
             by = "fgsea")
goplot(gse)

### write gse file
write.table(gse,file = "goGSEA.txt", sep= " ", 
            row.names = F)
goGSE <- read.table("goGSEA.txt")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + 
        facet_grid(.~.sign)

## KEGG pathway analysis

kk <- gseKEGG(genelist, 
        organism = "hsa",
        pvalueCutoff = 0.05,
        minGSSize = 3,
        maxGSSize = 800,
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")

dotplot(kk, showCategory=10, split=".sign") + 
        facet_grid(.~.sign)

head(kk)

##########################################################
## Again GSEA for CD133+FA
df1_CD133PFA$Log2FoldChange

gene.df_133PFA <- bitr(df1_CD133PFA$Gene.symbol,
                fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Hs.eg.db)

df1_CD133_ENSG <- rename(gene.df_133PFA, 
                          Gene.symbol = SYMBOL)

df_CD133_entrez_Gene <- inner_join(df1_CD133PFA,df1_CD133_ENSG, 
                                   by = "Gene.symbol")

##############################################################

gene_to_test <- df_CD133_entrez_Gene$ENSEMBL

go_results <- enrichGO(gene = gene_to_test,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")
go_biologicalprocess <- as.data.frame(go_results)

fit_go <- plot(barplot(go_results, showCategory = 20))

############################################################
df1_EpCAMPFA$LogFoldChange

gene.df_EpCAMPFA <- bitr(df1_EpCAMPFA$Gene.symbol,
     fromType = "SYMBOL",
     toType = c("ENSEMBL"),
     OrgDb = org.Hs.eg.db,
     drop = TRUE)
df1_EpCAMPFA_ENSG <- rename(gene.df_EpCAMPFA, 
                            Gene.symbol = SYMBOL)
## join
df_EpCAMPFA_ENSG_Gene <- inner_join(df1_EpCAMPFA, 
                                    df1_EpCAMPFA_ENSG, 
                                    by = "Gene.symbol")

## goGSEA
EpCAMPFA_Gene_to_test <- df_EpCAMPFA_ENSG_Gene$ENSEMBL

EpCAMPFA_Go_result <- enrichGO(gene = EpCAMPFA_Gene_to_test,
         OrgDb = "org.Hs.eg.db",
         keyType = "ENSEMBL",
         ont = "BP")

EpCAMPFA_go_BP <- as.data.frame(EpCAMPFA_Go_result)

EpCAMPFA_fit <- plot(barplot(EpCAMPFA_Go_result, showCategory = 10))

#########################################################################































