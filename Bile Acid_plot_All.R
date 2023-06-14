## SP FA bile acid genes
df_sp = read.csv("sp_genelist.csv", header = TRUE)
## read Hallmark gmt file.
gmt_gene_symbol <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")

bile_acid_genelist <- gmt_gene_symbol %>% filter(term == "HALLMARK_BILE_ACID_METABOLISM")

bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
bile_SPFA <- inner_join(bile_acid_genelist1,df_sp, by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change <- 2^bile_SPFA$Log2FoldChange

bile_SPFA1 <- cbind(bile_SPFA, FoldChange = Fold_Change)

## reorder the genes on based of fold change.
bile_gene <- bile_SPFA1[order(bile_SPFA1$FoldChange, decreasing = TRUE), ]

spfa_plot <- ggplot(bile_gene,aes(x = reorder(Gene.symbol,-FoldChange), 
                                  y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in SP FA")+
  theme_classic()
ggsave("spfa.jpg", plot = spfa_plot, width = 8, height = 6, dpi = 300)

######################################################################



df1_CD133PFA = read_excel("All_genelist.xlsx", sheet = 2)
## join genes from both data
#bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
df1_cd133_bilegene <- inner_join(bile_acid_genelist1,df1_CD133PFA, 
                                 by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change_cd133pFA <- 2^df1_cd133_bilegene$Log2FoldChange


df1_CD133PFA_bileFC <- cbind(df1_cd133_bilegene, 
                             FoldChange = Fold_Change_cd133pFA)

## reorder the genes on based of fold change.
df1_CD133PFA_bileFC1 <- df1_CD133PFA_bileFC[order(df1_CD133PFA_bileFC$FoldChange, decreasing = TRUE), ]

cd133pfaplot <- ggplot(df1_CD133PFA_bileFC1,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
     geom_bar(stat = "identity", position = "dodge", fill = "black") +
     xlab("Gene") +
     ylab("Fold Change") +
     ggtitle("Bile acid genes in CD133+FA")+
     theme_classic()

ggsave("cd133pfa.jpg", plot = cd133pfaplot, width = 8, height = 6, dpi = 300)

################################################################################

df1_EpCAMPFA = read_excel("All_genelist.xlsx", sheet = 3)
## join genes from both data
#bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
df1_EpCAMPFA_bilegene <- inner_join(bile_acid_genelist1,df1_EpCAMPFA, 
                                 by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change_EpCAMPFA <- 2^df1_EpCAMPFA_bilegene$LogFoldChange


df1_EpCAMPFA_bileFC <- cbind(df1_EpCAMPFA_bilegene, 
                             FoldChange = Fold_Change_EpCAMPFA)

## reorder the genes on based of fold change.
df1_EPCAMPFA_bileFC1 <- df1_EpCAMPFA_bileFC[order(df1_EpCAMPFA_bileFC$FoldChange, decreasing = TRUE), ]

# Split genes into up-regulated and down-regulated groups
# EpCAM_upregulated <- subset(df1_EPCAMPFA_bileFC1, FoldChange > 0)
# EpCAM_downregulated <- subset(df1_EPCAMPFA_bileFC1, FoldChange < 0)


plot_epcampfa_upreg <- ggplot(EpCAM_upregulated,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in EpCAM+FA")+
  theme_classic()


# plot_epcampfa_dnreg <- ggplot(EpCAM_downregulated,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
#   geom_bar(stat = "identity", position = "dodge", fill = "black") +
#   xlab("Gene") +
#   ylab("Fold Change") +
#   ggtitle("Bile acid genes in EpCAM+FA")+
#   theme_classic()
# Arrange the plots side by side
Epcampfa_plot <- plot_epcampfa_upreg + plot_epcampfa_dnreg


ggsave("epcampfa.jpg", plot = plot_epcampfa_upreg, width = 12, height = 6, dpi = 300)

################################################################################

df1_NSPFA = read_excel("All_genelist.xlsx", sheet = 4)
## join genes from both data
#bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
df1_NSPFA_bilegene <- inner_join(bile_acid_genelist1,df1_NSPFA, 
                                 by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change_NSPFA <- 2^df1_NSPFA_bilegene$Log2FoldChange


df1_NSPPFA_bileFC <- cbind(df1_NSPFA_bilegene, 
                             FoldChange = Fold_Change_NSPFA)

## reorder the genes on based of fold change.
df1_NSPPFA_bileFC1 <- df1_NSPPFA_bileFC[order(df1_NSPPFA_bileFC$FoldChange, decreasing = TRUE), ]

NSPFAplot <- ggplot(df1_NSPPFA_bileFC1,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in NSPFA")+
  theme_classic()

ggsave("NSPFAplot.jpg", plot = NSPFAplot, width = 6, height = 6, dpi = 300)

################################################################################


df1_CD133NFA = read_excel("All_genelist.xlsx", sheet = 5)
## join genes from both data
#bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
df1_CD133NFA_bilegene <- inner_join(bile_acid_genelist1,df1_CD133NFA, 
                                 by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change_CD133NFA <- 2^df1_CD133NFA_bilegene$Log2FoldChange


df1_CD133NFA_bileFC <- cbind(df1_CD133NFA_bilegene, 
                           FoldChange = Fold_Change_CD133NFA)

## reorder the genes on based of fold change.
df1_CD133NFA_bileFC1 <- df1_CD133NFA_bileFC[order(df1_CD133NFA_bileFC$FoldChange, decreasing = TRUE), ]

CD133NFAplot <- ggplot(df1_CD133NFA_bileFC1,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in CD133-FA")+
  theme_classic()

ggsave("CD133NFA.jpg", plot = CD133NFAplot, width = 6, height = 6, dpi = 300)

###############################################################################


df1_EpCAMNFA = read_excel("All_genelist.xlsx", sheet = 6)
## join genes from both data
#bile_acid_genelist1 <- rename(bile_acid_genelist, Gene.symbol = gene)
df1_EpCAMNFA_bilegene <- inner_join(bile_acid_genelist1,df1_EpCAMNFA, 
                                    by = "Gene.symbol")

## calculate fold change and add the column in data frame.
Fold_Change_EpCAMNFA <- 2^df1_EpCAMNFA_bilegene$Log2FoldChange


df1_EpCAMNFA_bileFC <- cbind(df1_EpCAMNFA_bilegene , 
                             FoldChange = Fold_Change_EpCAMNFA)

## reorder the genes on based of fold change.
df1_EpCAMNFA_bileFC1 <- df1_EpCAMNFA_bileFC[order(df1_EpCAMNFA_bileFC$FoldChange,
                                                  decreasing = TRUE), ]

# Split genes into up-regulated and down-regulated groups
 EpCAM_NFA_upreg <- subset(df1_EpCAMNFA_bileFC1, FoldChange > 1)
 EpCAM_NFA_dnreg <- subset(df1_EpCAMNFA_bileFC1, FoldChange < 1)

 EpCAM_NFA_upreg_plot <- ggplot(EpCAM_NFA_upreg,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "green") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in EpCAM-FA")+
  theme_classic()
 
 EpCAM_NFA_dnreg_plot <- ggplot(EpCAM_NFA_dnreg,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
   geom_bar(stat = "identity", position = "dodge", fill = "red") +
   xlab("Gene") +
   ylab("Fold Change") +
   ggtitle("Bile acid genes in EpCAM-FA")+
   theme_classic()
 
 EpCAM_NFA_final <- EpCAM_NFA_upreg_plot + EpCAM_NFA_dnreg_plot

ggsave("EpCAMNFA.jpg", plot =  EpCAM_NFA_final, width = 6, height = 6, dpi = 300)
