## differential genes-----------------------
setwd("~/gastroenterology")
pbmc <- readRDS("result/merged.rds")
celltype_color <- readRDS("color/celltype_color.rds")

pbmc <- subset(pbmc, group!="SC")
table(pbmc$group)

Epi <- subset(pbmc, celltype=="Epithelial_cell"&location!="Stomach")
table(Epi$group)
table(Epi$location)

comparelist <- list(c("HC","CEAS"))

Idents(Epi) <- "group"
markers <- FindMarkers(Epi, ident.1 = "CEAS",ident.2 = "HC",test.use = "wilcox")
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05)
data$cluster[data$regulate=="Up"]<-"CEAS"
data$cluster[data$regulate=="Down"]<-"HC"
data <- data[!is.na(data$cluster),]
data$gene <- rownames(data)
write.table(data,file="result/markers_Epi.txt",quote=F,sep="\t",row.names=F,col.names=T)

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
HC_markers <- data[data$regulate=="Down",]
HC_markers <- HC_markers[order(HC_markers$log2FoldChange),]

pdf("plot_new/volcano_Epi.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL", 
          label = "row", output = FALSE, custom_label = c(head(CEAS_markers$row,10),head(HC_markers$row,10)))
dev.off()

write.table(data[data$regulate!="Normal",], file = "result/Epi_markers.tsv", quote = F, sep = "\t")

sce.markers <- data[data$regulate!="Normal",]
sce.markers$cluster <- sce.markers$regulate
sce.markers$cluster[sce.markers$cluster=="Up"] <- "CEAS"
sce.markers$cluster[sce.markers$cluster=="Down"] <- "HC"
sce.markers$gene <- rownames(sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

## KEGG
xx <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", 
                     pvalueCutoff=0.01,pAdjustMethod = "BH",qvalueCutoff = 0.05) 
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

## GO
xx <- compareCluster(gcSample,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.01, ont="BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
p <- dotplot(xx)
pdf("plot_new/GO_Epi.pdf")
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

Epi <- NormalizeData(Epi)
Epi <- FindVariableFeatures(Epi)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Epi <- CellCycleScoring(Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Epi <- ScaleData(Epi, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Epi))

library(harmony)
Epi$batch <- Epi$sample_ident
Epi <- Epi %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

Epi <- Epi %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=0.5)
table(Epi$seurat_clusters)
Epi <- Epi %>%
  RunUMAP(dims = 1:30, reduction = "harmony")
DimPlot(Epi, label = T)
DimPlot(Epi, label = T, group.by = "group")

Idents(Epi) <- "seurat_clusters"
markers <- FindAllMarkers(Epi, only.pos = T, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

genes_to_check = unique(c('CHGA', 'TPH1', 'CES1',  'SLC38A11', 'RAB3C', #Enterochromaffin cells
                   'BEST4', 'CA7', 'CA4', 'SPIB', 'OTOP2', 'NOTCH2NL', #Enterocytes BEST4+ 
                   'CA1', 'SLC26A2', 'CA2', 'SLC26A3', 'KRT19', 'SELENBP1', 'PKIB', 'UGT2B17', 'CES2', #Enterocytes CA1+ CA2+ CA4-
                   'TMIGD1', 'MEP1A', 'APOA4', 'APOC3', 'APOA1', 'FABP6', # Enterocytes TMIGD1+ MEP1A+
                   'GSTA1', 'GSTA2', 'TMIGD1', 'MEP1A', #Enterocytes TMIGD1+ MEP1A+ GSTA1+
                   'PCSK1N', 'PYY', 'CHGA', 'GCG', 'CRYBA2', 'SCGN', 'FEV', 'SCG5', 'INSL5', 'MS4A8', #Enteroendocrine cells
                   'HBB', 'HBA2', 'HBA1', #Epithelial cells HBB+ HBA+
                   'MAFB', 'METTL12', #Epithelial cells METTL12+ MAFB+
                   'UBE2C', 'PTTG1', 'HMGB2', 'TOP2A', 'CKS2', 'CENPW', 'CDKN3', 'STMN1', 'TUBB4B', 'HIST1H4C', #Epithelial Cycling cells
                   'MUC2', 'RETNLB', 'SPINK4', 'ITLN1', 'CLCA1', 'FCGBP', 'TFF3', 'ST6GALNAC1', 'LRRC26', 'REP15', #Goblet cells MUC2+ TFF1- 
                   'MUC2', 'SPINK4', 'FCGBP', 'CLCA1', 'ZG16', 'TFF1', 'BCAS1', 'CEACAM5', #Goblet cells MUC2+ TFF1+ 
                   'SPINK4', 'MUC2', 'FCGBP', 'CLCA1', 'ITLN1', 'TFF3', 'TFF1', 'S100P', 'RETNLB', 'LRRC26', #Goblet cells SPINK4+ 
                   'DEFA5', 'DEFA6', 'REG3A', 'PRSS1', 'ITLN2', 'PLA2G2A', #Paneth cells
                   'OLFM4', 'REG1A', #Stem cells OLFM4+
                   'FABP1', 'GSTA1', 'AKR1C3', 'KRT19', 'MAOA', 'CES2', 'CBR1', 'RBP2', 'PTGR1', 'LIMA1', #Stem cells OLFM4+ GSTA1+
                   'LGR5', 'OLFM4', #Stem cells OLFM4+ LGR5+
                   'PCNA', 'RANBP1', 'OLFM4', 'STRA13', 'DUT', 'SIVA1', #Stem cells OLFM4+ PCNA+ 
                   'SH2D6', 'LRMP', '7SK_ENSG00000260682', 'AVIL', 'BMX', 'AZGP1', 'MATK', 'TRPM5'#Tuft cells
))
pdf("plot_new/heatmap_Epi_clusters.pdf", height = 10, width = 6)
DotPlot(Epi, features = genes_to_check)+coord_flip() + #翻转
  theme(panel.grid = element_blank(),
        axis.text.y=element_text(size = 8))+ #轴标签
  labs(x=NULL,y=NULL) +
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()


markers_df <- markers[markers$p_val_adj<0.01,]
ids=bitr(markers_df$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
markers_df=merge(markers_df,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(markers_df$ENTREZID, markers_df$cluster) 
lengths(gcSample)

## GO
xx <- compareCluster(gcSample,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.01, ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
p <- dotplot(xx)
pdf("plot_new/GO_Epi_cluster.pdf", height = 10, width = 8)
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

Epi$subtype <- "Unknown"
Epi$subtype[Epi$seurat_clusters%in%c(1,2,4)] <- "Enterocytes TMIGD1+ MEP1A+"
Epi$subtype[Epi$seurat_clusters%in%c(0,10)] <- "Goblet cells MUC2+ TFF1-"
Epi$subtype[Epi$seurat_clusters%in%c(11)] <- "Goblet cells MUC2+ TFF1+"
#Epi$subtype[Epi$seurat_clusters%in%c(16)] <- "Enteroendocrine cells"
Epi$subtype[Epi$seurat_clusters%in%c(6)] <- "Enterocytes CA1+ CA2+ CA4-"
Epi$subtype[Epi$seurat_clusters%in%c(8)] <- "Tuft cells"
Epi$subtype[Epi$seurat_clusters%in%c(9)] <- "Paneth cells"
Epi$subtype[Epi$seurat_clusters%in%c(3,5)] <- "Stem cells OLFM4+ GSTA1+"
DimPlot(Epi,group.by = "subtype")
SubtypeColor <- ArchRPalettes$stallion[1:length(unique(Epi$subtype))]
names(SubtypeColor) <- unique(Epi$subtype)
SubtypeColor["Unknown"]<-"grey"

SampleColor <- readRDS("color/SampleColor.rds")
GroupColor <- readRDS("color/GroupColor.rds")
LocationColor <- readRDS("color/LocationColor.rds")
PatientColor <- readRDS("color/PatientColor.rds")
pdf(file="plot_new/Epi.pdf")
DimPlot(Epi, label = T)
DimPlot(Epi, group.by = "sample_ident", cols = SampleColor, shuffle = T)
DimPlot(Epi, group.by = "group", cols = GroupColor, shuffle = T)
DimPlot(Epi, group.by = "patient_ident", cols = PatientColor, shuffle = T)
DimPlot(Epi, group.by = "location", cols = LocationColor, shuffle = T)
DimPlot(Epi, group.by = "Phase", shuffle = T)
DimPlot(Epi, group.by = "subtype", cols = SubtypeColor, shuffle = T)
ggplot(as.data.frame(table(Epi$sample_ident)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(Epi$group)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = GroupColor)+theme_bw()+xlab("Group")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(Epi$location)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = LocationColor)+theme_bw()+xlab("Location")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
dev.off()

table(Epi$subtype)
library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- Epi@meta.data
P1=CellInfo %>% ggplot(aes(x=group, fill=seurat_clusters)) +
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(seurat_clusters , fill=seurat_clusters))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)
P3=CellInfo %>% ggplot(aes(x=group, fill=subtype)) +
  geom_bar(color="black",position = "fill",width = 0.7) +scale_fill_manual(values=SubtypeColor)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P3g=ggplotGrob(P3)
P4=ggplot(CellInfo, aes(subtype , fill=subtype))+geom_bar(stat="count",colour = "black",width = 0.7)+scale_fill_manual(values=SubtypeColor)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P4g=ggplotGrob(P4)
pdf("plot_new/proportion_cluster_Epi.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
grid.arrange(grobs=list(P3g,P4g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

features<-list(artery1=c("HEY1","CXCL12", "AMD1", "FN1"),
               artery2=c("HEY1","CLU", "ELN", "IGFBP3"),
               capillary1=c("KDR", "CD36", "ID2", "SPP1"),
               capillary2=c("KDR", "CD36", "CA4", "FABP4", "LPL", "CLDN5", "FABP5"),
               capillary_venous=c("ACKR1", "CD36", "ICAM1", "CCL2", "SELE"),
               vein1=c("ACKR1", "RGCC", "POSTN", "VCAN", "NR2F2", "C7", "CYP1B1", "NR2F2", "APLNR"),
               vein2=c("ACKR1", "LAMA5"),
               vein3=c("ACKR1", "HLA-DQA2", "HLA-DQA1", "VCAM1", "PLAT"),
               vein4=c("ACKR1", "EDN1", "CYTL1"),
               angiogenic=c("KDR", "RGCC", "EDNRB", "VWA1", "ESM1", "INSR", "COL4A2", "COL4A1", "MCAM", "INSR", "NOTCH4"),
               lymphatic=c("CCL21", "PROX1", "LYVE", "PDPN", "TFF3"))

