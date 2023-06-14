library(msigdbr)
library(enrich)
all_gene_sets <-  msigdbr(species = "Homo sapiens")
head(all_gene_sets)

h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
head(h_gene_sets)

all_gene_sets %>%
  dplyr::filter(gs_cat == "H") %>%
  head()

msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
enricher(gene = genelist1, TERM2GENE = msigdbr_t2g)

gmt_gene_symbol <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")

bile_acid_genelist <- gmt_gene_symbol %>% filter(term == "HALLMARK_BILE_ACID_METABOLISM")

#gene.df <- rename(gene.df,Gene.symbol= SYMBOL)
#df_sp_entrez <- inner_join(df_sp,gene.df, by = "Gene.symbol")

bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
bile_SPFA <- inner_join(bile_acid_genelist1,df_sp, by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change <- 2^bile_SPFA$Log2FoldChange

bile_SPFA1 <- cbind(bile_SPFA, FoldChange = Fold_Change)

## reorder the genes on based of fold change.
bile_gene <- bile_SPFA1[order(bile_SPFA1$FoldChange, decreasing = TRUE), ]



















