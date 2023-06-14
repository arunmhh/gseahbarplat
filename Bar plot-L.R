rm(list = ls())
setwd("~/Desktop/Rstudio/Bar plot-enrichment analysis")

#---------------BAR PLOT-----------------
library(tidyverse)
library(ggplot2)
library(patchwork)

dat1 <- read.csv("Adeno_Top100_up_Pathway.csv", header = T)
dat1 <- dat1 %>% arrange(desc(LogP))
dat1$Term <- factor(dat1$Term, levels = rev(dat1$Term)) #🔧FIX🔧#

dat2 <- read.csv("Adeno_Top100_down_Pathway.csv", header = T)
dat2 <- dat2 %>% arrange(desc(LogP))
dat2$Term <- factor(dat2$Term, levels = rev(dat2$Term)) #🔧FIX🔧#

dat3 <- read.csv("Squamous_Top100_up_Pathway.csv", header = T)
dat3 <- dat3 %>% arrange(desc(LogP))
dat3$Term <- factor(dat3$Term, levels = rev(dat3$Term)) #🔧FIX🔧#

dat4 <- read.csv("Squamous_Top100_down_Pathway.csv", header = T)
dat4 <- dat4 %>% arrange(desc(LogP))
dat4$Term <- factor(dat4$Term, levels = rev(dat4$Term)) #🔧FIX🔧#

p1 <- ggplot(dat1, aes(y = Count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Adenocarcinoma Top 100 upregulated genes")+   #记得改名！！！！！！！！！！
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), #加在element_text()里面：标题位置：hjust = 0.5（中间）；1（右对齐）,字体：family = "Arial"；加粗：face = "bold"
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p2 <- ggplot(dat2, aes(y = Count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Adenocarcinoma Top 100 downregulated genes")+   #记得改名！！！！！！！！！！
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), #加在element_text()里面：标题位置：hjust = 0.5（中间）；1（右对齐）,字体：family = "Arial"；加粗：face = "bold"
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p3 <- ggplot(dat3, aes(y = Count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Squamous cell carcinoma Top 100 upregulated genes")+   #记得改名！！！！！！！！！！
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), #加在element_text()里面：标题位置：hjust = 0.5（中间）；1（右对齐）,字体：family = "Arial"；加粗：face = "bold"
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p4 <- ggplot(dat4, aes(y = Count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Squamous cell carcinoma Top 100 downregulated genes")+   #记得改名！！！！！！！！！！
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), #加在element_text()里面：标题位置：hjust = 0.5（中间）；1（右对齐）,字体：family = "Arial"；加粗：face = "bold"
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 2)


#ggsave(p,filename = "GO.pdf",width = 10,height = 7,dpi=300)
#ggsave(p,filename = "GO.jpg",width = 10,height = 7,dpi=300)

