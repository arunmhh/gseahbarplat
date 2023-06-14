library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggvenn)
library(qpcR)
## load data from excel table 
##(Header of table must be name of condition)
listsp <- read_excel("venn_gene_list.xlsx", sheet = 1)
listcd133 <- read_excel("venn_gene_list.xlsx", sheet = 2)
listepcam <- read_excel("venn_gene_list.xlsx", sheet = 3)
####################################################################
## convert different length of gene list into data frame using qpcR library

data_gene_venn <- qpcR:::cbind.na(listsp, listcd133,listepcam)

# Convert the data frame to a list
gene_df_list <- as.list(data_gene_venn)

## Visulize the venn diagram
upreg_venn <- ggvenn(gene_df_list, 
       fill_color = c("blue","yellow","green"),
       stroke_size = 1,
       show_percentage = FALSE,
       set_name_color = "black")

upreg_venn1 <- upreg_venn+
  labs(title = "Venn comparison of UP Regulated DEGs")+
  theme(plot.title = element_text(color="black", face="bold", size = 20,
                                  hjust= 0.5))

ggsave("venn_comparsion_msigdb_onco_SPCD133Epcam/venn_UpREg.jpg",
       plot = last_plot(),width = 10, height = 15, dpi = 300)

#####################################################################

raw_data <- read_excel("venn_gene_list.xlsx", sheet = 7)
colnames(raw_data) <- tolower(colnames(raw_data))

























