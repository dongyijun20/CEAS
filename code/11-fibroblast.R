## differential genes-----------------------
setwd("~/gastroenterology")
pbmc <- readRDS("result/merged.rds")
celltype_color <- readRDS("color/celltype_color.rds")

pbmc <- subset(pbmc, group!="SC")
table(pbmc$group)

fibro <- subset(pbmc, celltype=="Fibroblast")
table(fibro$group)
table(fibro$location)

comparelist <- list(c("HC","CEAS"))

Idents(fibro) <- "group"
markers <- FindMarkers(fibro, ident.1 = "CEAS",ident.2 = "HC",test.use = "wilcox")
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05)
data$cluster[data$regulate=="Up"]<-"CEAS"
data$cluster[data$regulate=="Down"]<-"HC"
data <- data[!is.na(data$cluster),]
data$gene <- rownames(data)
write.table(data,file="result/markers_fibro.txt",quote=F,sep="\t",row.names=F,col.names=T)

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
HC_markers <- data[data$regulate=="Down",]
HC_markers <- HC_markers[order(HC_markers$log2FoldChange),]

pdf("plot_new/volcano_fibro.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL", 
          label = "row", output = FALSE, custom_label = c(head(CEAS_markers$row,10),head(HC_markers$row,10)))
dev.off()

write.table(data[data$regulate!="Normal",], file = "result/fibro_markers.tsv", quote = F, sep = "\t")

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
pdf("plot_new/GO_fibro.pdf")
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

fibro <- NormalizeData(fibro)
fibro <- FindVariableFeatures(fibro)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fibro <- CellCycleScoring(fibro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# fibro <- ScaleData(fibro, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(fibro))

library(harmony)
fibro$batch <- fibro$sample_ident
fibro <- fibro %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

fibro <- fibro %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=0.5)
table(fibro$seurat_clusters)
fibro <- fibro %>%
  RunUMAP(dims = 1:30, reduction = "harmony")
DimPlot(fibro, label = T)
DimPlot(fibro, label = T, group.by = "group")

Idents(fibro) <- "seurat_clusters"
markers <- FindAllMarkers(fibro, only.pos = T, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

genes_to_check = unique(c('CCL11', 'ADAMDEC1', 'CCL13', 'HAPLN1', #Fibroblasts ADAMDEC1+
                          'KCNN3', 'LY6H', 'DPT', 'C7',  'SCN7A', #Fibroblasts KCNN3+ LY6H+ 
                          'NPY', 'SLITRK6', 'F3', 'EDNRB', 'NSG1', #Fibroblasts NPY+ SLITRK6+
                          'SLPI', 'SFRP2', 'IGFBP6', 'MFAP5', #Fibroblasts SFRP2+ SLPI+ 
                          'SMOC2', 'PTGIS', 'F3', 'PCSK6', 'ADAMTSL3', 'PCSK6'#Fibroblasts SMOC2+ PTGIS+
))
pdf("plot_new/heatmap_fibro_clusters.pdf", height = 6, width = 6)
DotPlot(fibro, features = genes_to_check)+coord_flip() + #翻转
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
pdf("plot_new/GO_fibro_cluster.pdf", height = 10, width = 8)
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

fibro$subtype <- "Fibroblasts ADAMDEC1+"
fibro$subtype[fibro$seurat_clusters%in%c(2)] <- "Fibroblasts NPY+ SLITRK6+"
DimPlot(fibro,group.by = "subtype")
SubtypeColor <- ArchRPalettes$stallion[1:length(unique(fibro$subtype))]
names(SubtypeColor) <- unique(fibro$subtype)
SubtypeColor["Unknown"]<-"grey"

SampleColor <- readRDS("color/SampleColor.rds")
GroupColor <- readRDS("color/GroupColor.rds")
LocationColor <- readRDS("color/LocationColor.rds")
PatientColor <- readRDS("color/PatientColor.rds")
pdf(file="plot_new/fibro.pdf")
DimPlot(fibro, label = T)
DimPlot(fibro, group.by = "sample_ident", cols = SampleColor, shuffle = T)
DimPlot(fibro, group.by = "group", cols = GroupColor, shuffle = T)
DimPlot(fibro, group.by = "patient_ident", cols = PatientColor, shuffle = T)
DimPlot(fibro, group.by = "location", cols = LocationColor, shuffle = T)
DimPlot(fibro, group.by = "Phase", shuffle = T)
DimPlot(fibro, group.by = "subtype", cols = SubtypeColor, shuffle = T)
ggplot(as.data.frame(table(fibro$sample_ident)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(fibro$group)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = GroupColor)+theme_bw()+xlab("Group")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(fibro$location)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = LocationColor)+theme_bw()+xlab("Location")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
dev.off()

table(fibro$subtype)
library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- fibro@meta.data
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
pdf("plot_new/proportion_cluster_fibro.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
grid.arrange(grobs=list(P3g,P4g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

