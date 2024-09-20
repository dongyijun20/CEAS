## differential genes-----------------------
setwd("~/gastroenterology")
pbmc <- readRDS("result/merged.rds")
celltype_color <- readRDS("color/celltype_color.rds")

pbmc <- subset(merge, group!="SC")
table(pbmc$group)

endo <- subset(pbmc, celltype=="Endothelial_cell")
table(endo$group)

comparelist <- list(c("HC","CEAS"))

Idents(endo) <- "group"
markers <- FindMarkers(endo, ident.1 = "CEAS",ident.2 = "HC",test.use = "wilcox")
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05)
data$cluster[data$regulate=="Up"]<-"CEAS"
data$cluster[data$regulate=="Down"]<-"HC"
data <- data[!is.na(data$cluster),]
data$gene <- rownames(data)
write.table(data,file="result/markers_endo.txt",quote=F,sep="\t",row.names=F,col.names=T)

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
HC_markers <- data[data$regulate=="Down",]
HC_markers <- HC_markers[order(HC_markers$log2FoldChange),]

pdf("plot_new/volcano_endo.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL",
          label = "row", output = FALSE, custom_label = c(CEAS_markers$row[1:10],HC_markers$row[1:10]))
dev.off()

write.table(data[data$regulate!="Normal",], file = "result/endo_markers.tsv", quote = F, sep = "\t")

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
pdf("plot_new/GO_endo.pdf")
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

endo <- NormalizeData(endo)
endo <- FindVariableFeatures(endo)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
endo <- CellCycleScoring(endo, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
endo <- ScaleData(endo, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(endo))

endo$batch <- endo$sample_ident
endo <- endo %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

endo <- endo %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=0.1)
table(endo$seurat_clusters)
endo <- endo %>%
  RunUMAP(dims = 1:30, reduction = "harmony")
DimPlot(endo, label = T)

Idents(endo) <- "seurat_clusters"
markers <- FindAllMarkers(endo, test.use = "wilcox")
markers_df = markers %>% group_by(cluster)# %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.table(markers_df,"result/markers_posandneg.tsv",quote = F,row.names = F)

pdf("plot_new/heatmap_endo_clusters.pdf")
DoHeatmap(subset(endo, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
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
pdf("plot_new/GO_endo_cluster.pdf", height = 10, width = 8)
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

df<-unique(ids$SYMBOL[match(unlist(str_split(xx@compareClusterResult$geneID
                                             [grep("vasculogenesis",xx@compareClusterResult$Description)],"/")),ids$ENTREZID)])
write.table(df,file="result/vasculogenesis.txt",quote=F,sep="\t",row.names=F, col.names = F)


SampleColor <- readRDS("color/SampleColor.rds")
GroupColor <- readRDS("color/GroupColor.rds")
LocationColor <- readRDS("color/LocationColor.rds")
PatientColor <- readRDS("color/PatientColor.rds")
pdf(file="plot_new/endo.pdf")
DimPlot(endo, label = T)
DimPlot(endo, group.by = "sample_ident", cols = SampleColor, shuffle = T)
DimPlot(endo, group.by = "group", cols = GroupColor, shuffle = T)
DimPlot(endo, group.by = "patient_ident", cols = PatientColor, shuffle = T)
DimPlot(endo, group.by = "location", cols = LocationColor, shuffle = T)
DimPlot(endo, group.by = "Phase", shuffle = T)
ggplot(as.data.frame(table(endo$sample_ident)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(endo$group)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = GroupColor)+theme_bw()+xlab("Group")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(endo$location)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = LocationColor)+theme_bw()+xlab("Location")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
dev.off()

## endothelial reference construction-----------------
seu.data <- read.csv("reference/GSE155109_bcTax_EC_enriched_raw_matrix.csv.gz")
rownames(seu.data) <- seu.data$Feature
seu.data <- seu.data[,-1]
seu <- CreateSeuratObject(counts = seu.data, project = "reference", min.cells = 1, min.features = 1)
meta <- read.csv("reference/GSE155109_bcTax_EC_enriched_metadata.csv.gz", sep = ";")
rownames(meta) <- meta$X
meta <- meta[,-1]
table(meta$cell_type_annotation)
seu <- AddMetaData(seu, as.data.frame(meta))
head(seu@meta.data)

# pre-process dataset (without integration)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:30, return.model = TRUE)

# map to reference
anchors <- FindTransferAnchors(reference = seu, query = endo, dims = 1:30,
                                        reference.reduction = "pca")
endo <- MapQuery(anchorset = anchors, reference = seu, query = endo,
                           refdata = list(EC_type = "cell_type_annotation"), 
                 reference.reduction = "pca", reduction.model = "umap")

#merge reference and query
seu$id <- 'reference'
endo$id <- 'query'
refquery <- merge(seu, endo)
refquery[["umap"]] <- merge(seu[["umap"]], endo[["ref.umap"]])

EC_color <- ArchRPalettes$kelly[1:length(unique(seu$cell_type_annotation))]
names(EC_color) <- unique(seu$cell_type_annotation)
saveRDS(EC_color, file = "color/EC_color.rds")

pdf("plot_new/ref.pdf", width = 5, height = 4)
DimPlot(seu, group.by = "cell_type_annotation", cols = EC_color)
DimPlot(endo, reduction = "umap", group.by = "predicted.EC_type", shuffle = TRUE, cols = EC_color) 
DimPlot(endo, reduction = "ref.umap", group.by = "predicted.EC_type", shuffle = TRUE, cols = EC_color)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)
# FeaturePlot(endo, features = c("prediction.score.lymphatic", "prediction.score.angiogenic", "prediction.score.activated_PCV"),  
#             reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
dev.off()

table(endo$predicted.EC_type)
library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- endo@meta.data
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
P3=CellInfo %>% ggplot(aes(x=group, fill=EC_type)) +
  geom_bar(color="black",position = "fill",width = 0.7) +scale_fill_nejm()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P3g=ggplotGrob(P3)
P4=ggplot(CellInfo, aes(EC_type , fill=EC_type))+geom_bar(stat="count",colour = "black",width = 0.7)+scale_fill_nejm()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P4g=ggplotGrob(P4)
pdf("plot_new/proportion_cluster_endo.pdf",width=7,height=5)
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

# ## ssGSEA
# exp=AverageExpression(endo,assays = "RNA")
# counts2=exp[["RNA"]]
# GSVA_hall <- gsva(expr=counts2,
#                   gset.idx.list=features,
#                   mx.diff=T, # 数据为正态分布则T，双峰则F
#                   kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
#                   parallel.sz=1) # 并行线程数目
# 
# head(GSVA_hall)
# 
# subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
# endo$EC_type <- sapply(endo$seurat_clusters,function(x) subtype[x])
# pheatmap::pheatmap(GSVA_hall, #热图的数据
#                    cluster_rows = T,#行聚类
#                    cluster_cols =T,#列聚类，可以看出样本之间的区分度
#                    show_colnames=T,
#                    scale = "column") #以行来标准化，这个功能很不错

endo$EC_type<-NA
endo$EC_type[endo$seurat_clusters==0|endo$seurat_clusters==2|endo$seurat_clusters==4]<-"capillary"
endo$EC_type[endo$seurat_clusters==1|endo$seurat_clusters==5|endo$seurat_clusters==6]<-"vein"
endo$EC_type[endo$seurat_clusters==3]<-"CEAS-related"
endo$EC_type[endo$seurat_clusters==7]<-"lymphatic"

pdf("plot_new/EC_cluster.pdf")
VlnPlot(endo, features = c("HEY1","KDR","CD36","ACKR1","CCL21","PROX1"), group.by = "seurat_clusters", stack = TRUE)
FeaturePlot(endo, features = c("HEY1","KDR","CD36","ACKR1","CCL21","PROX1"), max.cutoff = "q90", min.cutoff = "q10")
DimPlot(endo, group.by = "EC_type", shuffle = T)+scale_color_nejm()
dev.off()

## differential genes by celltype-----------------
Idents(endo) <- "EC_type"
markers <- FindAllMarkers(endo, only.pos = T, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("plot_new/heatmap_endo_types.pdf")
DoHeatmap(subset(endo, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
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
pdf("plot_new/GO_endo_types.pdf", height = 8, width = 9)
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

saveRDS(endo,"result/endothelial_cells_rmbatch.rds")

## expression of SLCO2A1 in EC types-------------
plot_df <- data.frame(expression=endo@assays$RNA@data["SLCO2A1",],EC_type=endo$EC_type, group=endo$group)
stat.test <- plot_df %>%
  group_by(EC_type) %>%
  t_test(expression ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(x = "EC_type", dodge = 0.8)
plot_df$EC_type<-as.factor(plot_df$EC_type)
bxp<-ggboxplot(
  plot_df, x = "EC_type", y = "expression", 
  color = "group", palette = GroupColor) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pdf("plot_new/SLCO2A1_expression_endo.pdf", height = 4, width = 6)
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0
  )+ggtitle("SLCO2A1 Expression in Endothelial Cells")
dev.off()

cap <- subset(endo, EC_type=="capillary1")
Idents(cap) <- "group"
markers <- FindAllMarkers(cap, only.pos = T, test.use = "wilcox")
markers_df <- markers[markers$p_val_adj<0.05,]
table(markers_df$cluster)

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
pdf("plot_new/GO_cluster_endo_cap.pdf", height = 6, width = 6)
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()




### co-expression
subset1 <- subset(endo, group=="CEAS")
subset1 <- subset1[,sample(colnames(subset1),1000)]
subset2 <- subset(endo, group=="HC")
subset2 <- subset2[,sample(colnames(subset2),1000)]
subset <- merge(subset1, subset2)

expr <- subset@assays$RNA@data
target <- expr[rownames(expr)=="SLCO2A1",]
expr <- expr[rownames(expr)!="SLCO2A1",]
cor_list <- sapply(1:nrow(expr), function(x) c(cor(expr[x,1:ncol(subset1)],target[1:ncol(subset1)]),
                                               cor(expr[x,(ncol(subset1)+1):(ncol(subset1)+ncol(subset2))],target[(ncol(subset1)+1):(ncol(subset1)+ncol(subset2))])))

df <- data.frame(gene = rownames(expr), cor_CEAS = cor_list[1,], cor_HC=cor_list[2,], 
                 cor_abnormal=cor_list[1,]/cor_list[2,])
df <- na.omit(df)
head(df)

write.table(df, file="result/correlation_outliers_endo.txt", quote = F, row.names = F)
df$cor_CEAS <- as.numeric(df$cor_CEAS)
df$cor_HC <- as.numeric(df$cor_HC)
select <- df$gene[abs(df$cor_CEAS)>abs(df$cor_HC)+0.1&df$cor_CEAS>0]
select

pdf("plot/correlation.pdf")
ggplot(data = df, aes(x=cor_HC, y=cor_CEAS)) +
  geom_point(size=0.3, color = ifelse(abs(df$cor_CEAS)>abs(df$cor_HC)+0.1&df$cor_CEAS>0, "red", "black"))+
  theme_classic(base_size = 16)+
  geom_text_repel(data=subset(df, abs(df$cor_CEAS)>abs(df$cor_HC)+0.1&df$cor_CEAS>0), aes(x=cor_HC, y=cor_CEAS, label=gene), 
                  min.segment.length = 0.05,segment.alpha=0.6, label.padding = 0.4, 
                  size=3,segment.colour = "grey50",
                  max.overlaps =30,nudge_x = -0.05,nudge_y=0.05)
dev.off()

wp2gene1 <- read.gmt("../neurosurgery/NES/reference/c2.cp.kegg.v7.5.1.symbols.gmt")
wp2gene2 <- read.gmt("../neurosurgery/NES/reference/h.all.v2023.1.Hs.symbols.gmt")

ewp1 <- enricher(select, TERM2GENE = wp2gene1, pvalueCutoff=0.05)
ewp2 <- enricher(select, TERM2GENE = wp2gene2, pvalueCutoff=0.05)

pdf("plot/outliers_function.pdf", height = 4, width = 7)
barplot(ewp1, showCategory=30, font.size = 10)
barplot(ewp2, showCategory=30, font.size = 10)
dev.off()

PPAR <- str_split(ewp1@result[grep("PPAR",ewp1@result$Description),"geneID"],"/")[[1]]
ewp1 <- enricher(PPAR, TERM2GENE = wp2gene1, pvalueCutoff=0.05)
ewp2 <- enricher(PPAR, TERM2GENE = wp2gene2, pvalueCutoff=0.05)

pdf("plot/PPAR_function.pdf", height = 4, width = 7)
barplot(ewp1, showCategory=30, font.size = 10)
barplot(ewp2, showCategory=30, font.size = 10)
dev.off()

Idents(pbmc) <- "celltype"
pdf("plot/PPAR_expression.pdf")
VlnPlot(pbmc, PPAR, pt.size = 0)
dev.off()

## 把outlier基因和差异基因取交集
inter <- intersect(select, sce.markers$gene)
ewp1 <- enricher(inter, TERM2GENE = wp2gene1, pvalueCutoff=1)
ewp2 <- enricher(inter, TERM2GENE = wp2gene2, pvalueCutoff=0.05)

pdf("plot/PPAR_function.pdf", height = 4, width = 7)
barplot(ewp1, showCategory=30, font.size = 10)
barplot(ewp2, showCategory=30, font.size = 10)
dev.off()

## ssGSEA
gene_sets <- list(wp2gene1$gene[grep("PPAR",wp2gene1$term)])
GSVA_hall <- gsva(expr=as.matrix(pbmc@assays$RNA@data),
                  gset.idx.list=gene_sets,
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目

head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=F,
                   scale = "column") #以行来标准化，这个功能很不错

### split by celltype
obj_list <- SplitObject(pbmc, split.by = "celltype")
markers_list <- lapply(obj_list, function(x){
  Idents(x) <- "group"
  FindMarkers(x, ident.1 = "CEAS",ident.2 = "HC",test.use = "wilcox")
})
markers_list_endo_up <- rownames(markers_list$Endothelial_cell)[markers_list$Endothelial_cell$p_val_adj<0.05]

## trajectory
library(Seurat)
library(patchwork)
library(monocle)

## monocle2
DefaultAssay(endo) <- 'RNA'
cds <- as.CellDataSet(endo)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
## ordering by marker gene per cluster
Idents(endo) <- "seurat_clusters"
deg <- FindAllMarkers(endo)
dim(deg)
deg = deg %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
sel.gene <- deg$gene

# cds <- detectGenes(cds, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(cds),
#                                     num_cells_expressed >= 10))
# marker_diff <- markerDiffTable(cds[expressed_genes,],
#                                cth,
#                                residualModelFormulaStr = "~Media + num_genes_expressed",
#                                cores = 1)
# candidate_clustering_genes <-
#   row.names(subset(marker_diff, qval < 0.01))
# marker_spec <-
#   calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
# head(selectTopMarkers(marker_spec, 3))
# semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
# HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
# plot_ordering_genes(HSMM)

cds <- monocle::setOrderingFilter(cds, sel.gene)

## dimension reduciton
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

## ordering cells
cds <- monocle::orderCells(cds)

saveRDS(cds, file = "result/endo_cds.rds")

library(ggpubr)
df<-pData(cds)

pdf("plot_new/monocle2_endo.pdf")
monocle::plot_cell_trajectory(cds, color_by = "group", cell_size = 1) + scale_color_manual(values=GroupColor)
monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 1)
monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 0.5)  + facet_wrap(~seurat_clusters)
monocle::plot_cell_trajectory(cds, color_by = "State", cell_size = 1) + scale_color_jama()
monocle::plot_cell_trajectory(cds, color_by = "EC_type", cell_size = 0.5)+scale_color_nejm()
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1)
monocle::plot_complex_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 1)
monocle::plot_complex_cell_trajectory(cds, color_by = "State", cell_size = 1) + scale_color_nejm()
monocle::plot_complex_cell_trajectory(cds, color_by = "group", cell_size = 1) + scale_color_manual(values=GroupColor)
ggplot(df,aes(Pseudotime,color=seurat_clusters,fill=seurat_clusters))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()
ggplot(df,aes(Pseudotime,color=group,fill=group))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_manual(values=GroupColor)+scale_color_manual(values=GroupColor)
ggplot(df,aes(Pseudotime,color=State,fill=State))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_nejm()+scale_color_jama()
ggplot(df,aes(Pseudotime,color=EC_type,fill=EC_type))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_nejm()+scale_color_nejm()
dev.off()

pdf("plot_new/proportion_state_endo.pdf",width=6,height=5)
df %>% ggplot(aes(x=group, fill=State)) + scale_fill_nejm()+
  geom_bar(color="black",position = "fill",width = 0.7) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
df %>% ggplot(aes(x=State, fill=seurat_clusters)) +
  geom_bar(color="black",position = "fill",width = 0.7) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
df %>% ggplot(aes(x=State, fill=group)) +
  geom_bar(color="black",position = "fill",width = 0.7) + scale_fill_manual(values = GroupColor) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
dev.off()

plot_genes_in_pseudotime(cds[c("HEY1", "CXCL12", 
                               "ACKR1", "VCAM1", 
                               "CD36", "KDR", "EDNRB"),], ncol = 2, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res$gene_short_name[diff_test_res$qval<0.01]
length(sig_gene_names)
hm<-plot_pseudotime_heatmap(cds[sig_gene_names,],
                            num_clusters = 4, 
                            cores = 1, return_heatmap = T)

df_row_cluster = data.frame(cluster = cutree(hm$tree_row, k = 4))
df_row_cluster$gene <- rownames(df_row_cluster)
dim(df_row_cluster)

ids=bitr(df_row_cluster$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
df_row_cluster=merge(df_row_cluster,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(df_row_cluster$ENTREZID, df_row_cluster$cluster) 
lengths(gcSample)

## GO
xx <- compareCluster(gcSample,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.01, ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
p <- dotplot(xx)
pdf("plot_new/GO_cluster_endo_monocle2.pdf", height = 6, width = 8)
hm
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()


BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-5)),],
                            branch_point = 1,
                            num_clusters = 5,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

