endo <- readRDS("result/endothelial_cells_rmbatch.rds")
DimPlot(endo)
Idents(endo)<-"seurat_clusters"

endo$label <- paste(endo$group,endo$seurat_clusters,sep = "_")
table(endo$label)

interest <- subset(endo, label%in%c("CEAS_0","CEAS_3","HC_0"))
Idents(interest) <- "label"
markers1 <- FindMarkers(interest, ident.1 = "CEAS_0", ident.2 = "HC_0") # fold change 算的是ident.1/ident.2
markers2 <- FindMarkers(interest, ident.1 = "CEAS_3", ident.2 = "CEAS_0")
head(markers1)
head(markers2)

data <- add_regulate(markers1, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
control_markers <- data[data$regulate=="Down",]
control_markers <- control_markers[order(control_markers$log2FoldChange),]
data1 <- data

data <- add_regulate(markers2, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
CEAS_markers <- data[data$regulate=="Up",]
CEAS_markers <- CEAS_markers[order(CEAS_markers$log2FoldChange, decreasing = T),]
control_markers <- data[data$regulate=="Down",]
control_markers <- control_markers[order(control_markers$log2FoldChange),]
data2 <- data

library(ggVolcano)


pdf("plot_new/volcano_endo_subluster.pdf")
ggvolcano(data1, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL",
          label = "row", output = FALSE, custom_label = c(CEAS_markers$row[1:10],control_markers$row[1:10])) + ggtitle("Cluster0: CEAS vs HC")
ggvolcano(data2, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL",
          label = "row", output = FALSE, custom_label = c(CEAS_markers$row[1:10],control_markers$row[1:10])) + ggtitle("CEAS: Cluster3 vs Cluster0")
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
sce.markers <- data1[data1$regulate!="Normal",]
sce.markers$cluster <- sce.markers$regulate
sce.markers$cluster[sce.markers$cluster=="Up"] <- "CEAS_0"
sce.markers$cluster[sce.markers$cluster=="Down"] <- "HC_0"
sce.markers$gene <- rownames(sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01
)
p1 <- dotplot(xx)

sce.markers <- data2[data2$regulate!="Normal",]
sce.markers$cluster <- sce.markers$regulate
sce.markers$cluster[sce.markers$cluster=="Up"] <- "CEAS_3"
sce.markers$cluster[sce.markers$cluster=="Down"] <- "CEAS_0"
sce.markers$gene <- rownames(sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01
)
p2 <- dotplot(xx)

pdf("plot_new/GO_endo_subcluster.pdf")
p1 + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2 + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))  + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

### Seurat workflow
interest <- NormalizeData(interest)
interest <- FindVariableFeatures(interest)
interest <- ScaleData(interest)
interest <- RunPCA(interest)
interest <- RunUMAP(interest, dims = 1:40)
interest <- FindNeighbors(interest, dims = 1:40)
interest <- FindClusters(interest, resolution = 0.5)
DimPlot(interest, group.by = "label", label = T)
DimPlot(interest)

saveRDS(interest, file = "result/endo_subset.rds")

## trajectory-monocle3
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(viridis)

immune.combined <- interest
myeloid.cds <- as.cell_data_set(immune.combined)
myeloid.cds<-preprocess_cds(myeloid.cds)

# myeloid.cds <- cluster_cells(cds = myeloid.cds, reduction_method = "UMAP")
# myeloid.cds <- learn_graph(myeloid.cds, use_partition = TRUE)
# plot_cells(myeloid.cds)

myeloid.cds <- learn_graph(myeloid.cds, use_partition = F, close_loop = FALSE,
                           learn_graph_control = list(minimal_branch_len=18,euclidean_distance_ratio=5))
myeloid.cds <- order_cells(myeloid.cds, reduction_method = "UMAP")
plot_cells(myeloid.cds)
#hsc <- colnames(immune.combined)[which(immune.combined$seurat_clusters=='0')]
# order cells
#myeloid.cds <- order_cells(myeloid.cds, reduction_method = "UMAP", root_cells = hsc)
#plot_cells(myeloid.cds)

saveRDS(myeloid.cds,file="result/monocle3_endo.rds")

myeloid.cds <- readRDS("result/monocle3_endo.rds")
immune.combined <- AddMetaData(
  object = immune.combined,
  metadata = myeloid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "monocle3"
)

pdf(file="plot_new/monocle3.pdf", height = 4, width = 4)
# plot trajectories colored by pseudotime
plot_cells(
  cds = myeloid.cds,
  color_cells_by = "label",
  show_trajectory_graph = TRUE,
  label_leaves = TRUE, 
  label_branch_points = TRUE)
plot_cells(myeloid.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
FeaturePlot(immune.combined, c("monocle3"), pt.size = 0.1) & scale_color_viridis()
dev.off()

cds_subset <- choose_cells(myeloid.cds)
plot_cells(cds_subset)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
