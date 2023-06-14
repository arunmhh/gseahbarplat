## SP FA MsigDB top pathways
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)

sp_data <- read.csv("SP FA_MSigdb.csv", header = T)
sp_data <- sp_data %>% arrange(desc(LogP))
sp_data$Term <- factor(sp_data$Term, levels = rev(sp_data$Term))

spfa<- ggplot(sp_data, aes(y = count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("SP FA Pathways")+   
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

ggsave(spfa,"spfa_msigdb.jpg", width = 12, height = 10)

##############################
nsp_data <- read_excel("SP FA_MSigdb.xlsx", sheet = 2)
nsp_data <- nsp_data %>% arrange(desc(LogP))
nsp_data$Term <- factor(nsp_data$Term, levels = rev(nsp_data$Term))

nspfa<- ggplot(nsp_data, aes(y = count, x = Term, fill= LogP)) + 
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("NSPFA Pathways")+   
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), 
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

ggsave("nspfa_msigdb.jpg",plot = last_plot(), width = 12, height = 10, dpi = 300)
######################################################################
cd133nfa_data <- read_excel("SP FA_MSigdb.xlsx", sheet = 4)
cd133nfa_data <- cd133nfa_data %>% arrange(desc(LogP))
cd133nfa_data$Term <- factor(cd133nfa_data$Term, levels = rev(cd133nfa_data$Term))

cd133nfa<- ggplot(cd133nfa_data, aes(y = count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("CD133- FA Pathways")+   
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), 
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

ggsave("cd133nfa_msigdb.jpg",plot = last_plot(), width = 12, height = 10, dpi = 300)
######################################################################

epcampfa_data <- read_excel("SP FA_MSigdb.xlsx", sheet = 5)
epcampfa_data <- epcampfa_data %>% arrange(desc(LogP))
epcampfa_data$Term <- factor(epcampfa_data$Term, 
                             levels = rev(epcampfa_data$Term))

epcampfa<- ggplot(epcampfa_data, aes(y = count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("EpCAM+ FA Pathways")+   
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), 
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

ggsave("epcampfa_msigdb.jpg",plot = last_plot(), width = 12, height = 10,
       dpi = 300)

##############################################################################

epcamnfa_data <- read_excel("SP FA_MSigdb.xlsx", sheet = 6)
epcamnfa_data <- epcamnfa_data %>% arrange(desc(LogP))
epcamnfa_data$Term <- factor(epcamnfa_data$Term, 
                             levels = rev(epcamnfa_data$Term))

epcamnfa<- ggplot(epcamnfa_data, aes(y = count, x = Term, fill = LogP)) + 
  geom_bar(stat = "identity", position = "dodge")+
  ggtitle("EpCAM- FA Pathways")+   
  facet_grid(Category~., scales = "free", space = "free") + 
  scale_fill_gradient(low="#F11712",high="#0099F7") +
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = (20)), 
        strip.text.y = element_text(size = 14),
        legend.position ="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

ggsave("epcamnfa_msigdb.jpg",plot = last_plot(), width = 12, height = 10,
       dpi = 300)































































