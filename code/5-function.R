setwd("~/gastroenterology")
pbmc <- readRDS("result/merged_rmbatch.rds")
celltype_color <- readRDS("color/celltype_color.rds")

pbmc <- subset(merge, group!="SC")

plot_df <- data.frame(expression=pbmc@assays$RNA@data["SLCO2A1",],group=pbmc$group,celltype=pbmc$celltype)
stat.test <- plot_df %>%
  group_by(celltype) %>%
  t_test(expression ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp<-ggboxplot(
  plot_df, x = "celltype", y = "expression", 
  color = "group", palette = GroupColor) +
  ggtitle("SLCO2A1 Expression in All Cells") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf("plot_new/SLCO2A1_expression_celltype.pdf", height = 4, width = 7)
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0
)
dev.off()

### SLCO2A1 expression proportion
pbmc$SLCO2A1_exp <- ifelse(test = pbmc@assays$RNA@counts["SLCO2A1",]>0,yes = "1",no = "0")
table(pbmc$SLCO2A1_exp)
SLCO2A1 <- subset(pbmc, SLCO2A1_exp=="1")

library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- SLCO2A1@meta.data
#CellInfo <- CellInfo[CellInfo$sample_ident%in%c("Sample7","Sample8"),]
P1=CellInfo %>% ggplot(aes(x=group, fill=fct_rev(celltype))) +scale_fill_manual(values = celltype_color)+
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

pro_CEAS<-table(SLCO2A1[,SLCO2A1$group=="CEAS"]$celltype)/table(pbmc[,pbmc$group=="CEAS"]$celltype)
a<-table(SLCO2A1[,SLCO2A1$group=="HC"]$celltype)
a<-c(a,"NK_cell"=0,"Neutrophil"=0,"Smooth_muscle_cell"=0)
a<-a[names(pro_CEAS)]
pro_HC<-a/table(pbmc[,pbmc$group=="HC"]$celltype)
df <- data.frame(value=c(pro_CEAS,pro_HC),group=c(rep("CEAS",12),rep("HC",12)),celltype=c(df$CEAS.Var1,df$CEAS.Var1))

pdf("plot_new/proportion_SLCO2A1.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
ggplot(df, aes(celltype, value, fill=group))+geom_col(position="dodge")+  
  scale_fill_manual(values = GroupColor)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


epi <- subset(pbmc, celltype=="Epithelial_cell")
endo <- subset(pbmc, celltype=="Endothelial_cell")

epi
endo
Idents(endo) <- "group"
Idents(epi) <- "group"

comparelist <- list(c("HC","CEAS"))

pdf("plot/SLCO2A1_expression.pdf", height = 4, width = 4)
plot_df <- data.frame(expression=endo@assays$RNA@data["SLCO2A1",],group=endo$group)
ggplot(plot_df, aes(group, expression, fill = group)) +
  geom_violin() +
  scale_fill_manual(values = GroupColor) +
  geom_signif(comparisons = comparelist,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+ theme_classic()+
  ggtitle("SLCO2A1 Expression in Endothelial Cells")
plot_df <- data.frame(expression=epi@assays$RNA@data["SLCO2A1",],group=epi$group)
ggplot(plot_df, aes(group, expression, fill = group)) +
  geom_violin() +
  scale_fill_manual(values = GroupColor) +
  geom_signif(comparisons = comparelist,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+ theme_classic()+
  ggtitle("SLCO2A1 Expression in Epithelial Cells")
dev.off()

Idents(endo) <- "group"
markers <- FindMarkers(endo, ident.1 = "CEAS",ident.2 = "HC",test.use = "wilcox")

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
control_markers <- data[data$regulate=="Down",]
control_markers <- control_markers[order(control_markers$log2FoldChange),]

pdf("plot/volcano_epi.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL",
          label = "row", output = FALSE, custom_label = c(CEAS_markers$row[1:10],control_markers$row[1:10]))
dev.off()

write.table(data[data$regulate!="Normal",], file = "result/endo_markers.tsv", quote = F, sep = "\t")

sce.markers <- data[data$regulate!="Normal",]
sce.markers$cluster <- sce.markers$regulate
sce.markers$cluster[sce.markers$cluster=="Up"] <- "CEAS"
sce.markers$cluster[sce.markers$cluster=="Down"] <- "control"
sce.markers$gene <- rownames(sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

# ## KEGG
# xx <- compareCluster(gcSample,
#                      fun = "enrichKEGG", 
#                      organism = "hsa", pvalueCutoff = 0.05
# )
# p <- dotplot(xx)
# p + theme(axis.text.x = element_text(
#   angle = 45,
#   vjust = 0.5, hjust = 0.5
# ))

## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01
)
p <- dotplot(xx)
pdf("plot/GO_epi.pdf")
p + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

table(endo$batch)
table(epi$batch)
endo$batch <- endo$group
endo <- endo %>% 
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

endo <- endo %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 
endo <- endo %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(endo, group.by = "Phase")
DimPlot(endo, group.by = "group")

saveRDS(endo,"result/endothelial_cells_rmbatch.rds")

markers <- FindAllMarkers(endo, only.pos = T, test.use = "wilcox")

markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
write.table(markers_df,file="result/markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

markers_df <- read.table("result/markers200.txt", header = T)
#markers_df <- markers_df[-grep("^MT",markers_df$gene),]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

### co-expression
subset1 <- subset(endo, group=="CEAS")
subset1 <- subset1[,sample(colnames(subset1),1000)]
subset2 <- subset(endo, group=="control")
subset2 <- subset2[,sample(colnames(subset2),1000)]
subset <- merge(subset1, subset2)

expr <- subset@assays$RNA@data
target <- expr[rownames(expr)=="SLCO2A1",]
expr <- expr[rownames(expr)!="SLCO2A1",]
cor_list <- sapply(1:nrow(expr), function(x) c(cor(expr[x,1:ncol(subset1)],target[1:ncol(subset1)]),
                                               cor(expr[x,(ncol(subset1)+1):(ncol(subset1)+ncol(subset2))],target[(ncol(subset1)+1):(ncol(subset1)+ncol(subset2))])))

df <- data.frame(gene = rownames(expr), cor_CEAS = cor_list[1,], cor_control=cor_list[2,], 
                 cor_abnormal=cor_list[1,]/cor_list[2,])
df <- na.omit(df)
head(df)

write.table(df, file="result/correlation_outliers_endo.txt", quote = F, row.names = F)
df$cor_CEAS <- as.numeric(df$cor_CEAS)
df$cor_control <- as.numeric(df$cor_control)
select <- df$gene[abs(df$cor_CEAS)>abs(df$cor_control)+0.1&df$cor_CEAS>0]
select

pdf("plot/correlation.pdf")
ggplot(data = df, aes(x=cor_control, y=cor_CEAS)) +
  geom_point(size=0.3, color = ifelse(abs(df$cor_CEAS)>abs(df$cor_control)+0.1&df$cor_CEAS>0, "red", "black"))+
  theme_classic(base_size = 16)+
  geom_text_repel(data=subset(df, abs(df$cor_CEAS)>abs(df$cor_control)+0.1&df$cor_CEAS>0), aes(x=cor_control, y=cor_CEAS, label=gene), 
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
  FindMarkers(x, ident.1 = "CEAS",ident.2 = "control",test.use = "wilcox")
})
markers_list_endo_up <- rownames(markers_list$Endothelial_cell)[markers_list$Endothelial_cell$p_val_adj<0.05]





