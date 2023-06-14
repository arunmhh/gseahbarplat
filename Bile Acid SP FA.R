## blie acid sp FA.
  spfa_plot <- ggplot(bile_gene,aes(x = reorder(Gene.symbol,-FoldChange), y = FoldChange)) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  xlab("Gene") +
  ylab("Fold Change") +
  ggtitle("Bile acid genes in SP FA")+
    theme_classic()