merge <- readRDS("result/merged.rds")

markers <- FindAllMarkers(merge, only.pos = T, test.use = "wilcox")

markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
write.table(markers_df,file="result/markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

Idents(merge)<-"seurat_clusters"
pdf("plot_new/heatmap.pdf")
DoHeatmap(subset(merge, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
dev.off()

markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("plot_new/dotplot.pdf", height = 10, width = 6)
DotPlot(merge, features = unique(markers_df$gene)) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 8))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

table(merge$celltype_hpca_main, merge$seurat_clusters)
table(merge$celltype_bped_main, merge$seurat_clusters)

merge$celltype <- "Unknown"
merge$celltype[merge$seurat_clusters==0|merge$seurat_clusters==2] <- "T_cell"
merge$celltype[merge$seurat_clusters==1|merge$seurat_clusters==4|merge$seurat_clusters==5|merge$seurat_clusters==13] <- "B_cell"
merge$celltype[merge$seurat_clusters==3|merge$seurat_clusters==6|merge$seurat_clusters==8|merge$seurat_clusters==12|merge$seurat_clusters==20|merge$seurat_clusters==21] <- "Epithelial_cell"
merge$celltype[merge$seurat_clusters==9] <- "Neutrophil"
merge$celltype[merge$seurat_clusters==14] <- "Monocyte/Macrophage"
merge$celltype[merge$seurat_clusters==10] <- "Endothelial_cell"
merge$celltype[merge$seurat_clusters==15] <- "Mast_cell"
merge$celltype[merge$seurat_clusters==18] <- "Neuron"
merge$celltype[merge$seurat_clusters==7|merge$seurat_clusters==16] <- "Fibroblast"
merge$celltype[merge$seurat_clusters==19] <- "Smooth_muscle_cell"
merge$celltype[merge$seurat_clusters==11] <- "NK_cell"
merge$celltype[merge$seurat_clusters==17] <- "Proliferating_cell"
table(merge$celltype)

celltype_color <- ArchRPalettes$stallion[1:length(unique(merge$celltype))]
names(celltype_color) <- c(unique(merge$celltype))
saveRDS(celltype_color, file = "color/celltype_color.rds")

saveRDS(merge,"result/merged.rds")

pdf("plot_new/annotated.pdf")
DimPlot(merge, label=T)
FeaturePlot(merge, "SLCO2A1", min.cutoff = "q10", max.cutoff = "q90",cols = c("lightgrey" ,"#DE1F1F"))
DimPlot(merge, group.by = "celltype", cols = celltype_color)
dev.off()

pdf("plot_new/annotated_split.pdf", width = 10, height = 4)
DimPlot(merge,split.by = "group",group.by = "celltype", cols = celltype_color)
dev.off()

library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- merge@meta.data
CellInfo <- CellInfo[CellInfo$sample_ident%in%c("Sample7","Sample8"),]
P1=CellInfo %>% ggplot(aes(x=sample_ident, fill=fct_rev(celltype))) +scale_fill_manual(values = celltype_color)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(celltype , fill=celltype))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  scale_fill_manual(values = celltype_color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)
pdf("plot_new/proportion_celltype_sample.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

epi <- subset(merge, celltype=="Epithelial_cell")

# view cell cycle scores and phase assignments
pdf("plot_new/cell_cycle.pdf")
DimPlot(merge, group.by = "Phase")
dev.off()

markers_df <- read.table("result/markers200.txt", header = T)
#markers_df <- markers_df[-grep("^MT",markers_df$gene),]
top20.markers = markers_df %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
set.seed(12)
CellInfo <- merge@meta.data
Idents(merge)=CellInfo$sample_ident
#sampled.cells <- sample(x = merge@active.ident, size = 5000, replace = F)
subset<-subset(merge, group!="SC")
subset=subset(subset, downsample=1000)
merge.sce=as.SingleCellExperiment(subset)
library(SingleCellExperiment)
plot.data<-as.data.frame(assay(merge.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
CellInfoS=subset@meta.data

column_annot <-CellInfoS[,c("celltype","group"),drop=F]
column_annot$group = as.factor(as.character(column_annot$group))
column_annot=with(column_annot, column_annot[order(group), , drop=F])
column_annot=with(column_annot, column_annot[order(celltype), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]

SampleColors<-readRDS("color/GroupColor.rds")
AssignmentColors<-readRDS("color/celltype_color.rds")
# library(ArchR)
# ClusterColors<-ArchRPalettes$stallion[1:17]
# names(ClusterColors)<-c(0:16)
column.colors=list()
column.colors[["group"]]<-SampleColors
column.colors[["celltype"]]<-AssignmentColors
sample=as.matrix(column_annot[,c("group"),drop=F])
seurat_clusters=as.matrix(column_annot[,c("celltype"),drop=F])

library(ComplexHeatmap)
colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
genes= top2$gene
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$group)),title="group",legend_gp = gpar(fill=GroupColor[levels(as.factor(column_annot$group))],fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$celltype)),title="celltype",legend_gp = gpar(fill=cluster_color[levels(as.factor(column_annot$celltype))],fontsize=5))
#lgd3=Legend(labels = levels(as.factor(column_annot$Fragment)),title="Fragment",legend_gp = gpar(fill=column.col3,fontsize=5))

pdf("plot_new/heatmap_top20_seurat_clusters.pdf",width=7.5,height=9)
draw(HM,heatmap_legend_list = list(lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()



