color_database <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
  '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
  '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
  '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
  '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
  '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
  '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
  '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
  '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
  '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
  '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
  '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
saveRDS(color_database, file = "color/color_database.rds")

## prepare query data, colon stromal cells--------
merge <- readRDS("result/merged.rds")
DimPlot(merge, group.by="celltype")
merge_str <- subset(merge, celltype %in% c("Endothelial_cell","Fibroblast","Neuron","Smooth_muscle _cell") &
                      location == "Colon")
merge_str <- RunUMAP(merge_str, dims = 1:30)
DimPlot(merge_str, group.by = "celltype")

## prepare reference data, scIBD--------
scIBD <- readRDS("reference/scIBD/scIBD.gex_matrix.rds")
meta_scIBD <- read.csv("reference/scIBD/scIBD_sample_meta_data.txt", sep = "\t")
meta_scIBD <- column_to_rownames(meta_scIBD, "cell_id")

scIBD_subset <- subset(scIBD, disease %ni% c("Colitis_inflamed","CD_inflamed/CD_non_inflamed") & tissue.sub %in% c("duodenum","ileum","colon") & stage == "adult")
meta_scIBD_subset <- meta_scIBD[colnames(scIBD_subset),]
umap_reduction <- CreateDimReducObject(embeddings = as.matrix(meta_scIBD_subset[,c("gUMAP_1","gUMAP_2")]), key = "gUMAP_", assay = DefaultAssay(scIBD))
scIBD_subset@reductions$umap <- umap_reduction

# set UMAP models
umap_new_model <- list()
umap_new_model$n_epochs <- 500
umap_new_model$alpha <-1
umap_new_model$method <- "umap"
umap_new_model$negative_sample_rate <- 5
umap_new_model$gamma <- 1
umap_new_model$approx_pow <- 0
umap_new_model$n_neighbors <- 30
umap_new_model$metric$cosine <- list()
umap_new_model$embedding <- scIBD_subset[["umap"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap_new_model$a <- ab_param["a"]
umap_new_model$b <- ab_param["b"]
scIBD_subset[["umap"]]@misc$model <- umap_new_model

# annotate clusters
scIBD_subset$cluster <- case_when(scIBD_subset$major_cluster %in% c("Myeloid","CD4T","CD8T","ILC","B_Plasma") ~ "Immune",
                                  scIBD_subset$major_cluster %in% c("Mesenchymal","Neural","Endothelial") ~ "Stromal",
                                  scIBD_subset$major_cluster %in% c("Epithelial") ~ "Epithelial")

scIBD_subset <- NormalizeData(scIBD_subset)
scIBD_subset <- FindVariableFeatures(scIBD_subset)
scIBD_subset <- ScaleData(scIBD_subset)
scIBD_subset <- RunPCA(scIBD_subset, npcs = 50)

DimPlot(scIBD_subset, group.by = "cluster")
DimPlot(scIBD_subset, group.by = "major_cluster")
DimPlot(scIBD_subset, group.by = "minor_cluster")

saveRDS(scIBD_subset, file = "reference/scIBD/scIBD_subset4use.rds")

# select two technologies for the query datasets
anchors <- FindTransferAnchors(reference = scIBD_subset, query = merge, dims = 1:30, 
                               reference.reduction = "pca")
merge <- MapQuery(
  anchorset = anchors,
  query = merge,
  reference = scIBD_subset,
  refdata = list(
    cluster = "cluster",
    major_cluster = "major_cluster",
    minor_cluster = "minor_cluster"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(merge, reduction = "ref.umap", group.by = "predicted.cluster", label = TRUE, label.size = 3)
DimPlot(merge, reduction = "ref.umap", group.by = "predicted.major_cluster", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(merge, reduction = "ref.umap", group.by = "predicted.minor_cluster")

DimPlot(merge, reduction = "umap", group.by = "predicted.cluster", label = TRUE, label.size = 3)
DimPlot(merge, reduction = "umap", group.by = "predicted.major_cluster", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(merge, reduction = "umap", group.by = "predicted.minor_cluster")

# correct by getting the most frequent celltype for each seurat cluster
merge <- FindClusters(merge, resolution = 10)
DimPlot(merge)

correct_cluster <- function(seurat_obj, correct_group){
  # Extract the metadata
  metadata <- seurat_obj@meta.data
  # Get the most frequent cell type for each cluster
  table_cluster <- table(seurat_obj$seurat_clusters, as.vector(seurat_obj[[correct_group]])[[1]])
  cluster_celltype <- apply(table_cluster,1,function(x) colnames(table_cluster)[which.max(x)])
  cluster_vector <- sapply(seurat_obj$seurat_clusters, function(x) cluster_celltype[x])
  names(cluster_vector) <- sub("(.*)\\..*","\\1",names(cluster_vector))
  return(cluster_vector)
}

merge$corrected.cluster <- correct_cluster(seurat_obj = merge, correct_group = "predicted.cluster")
merge$corrected.major_cluster <- correct_cluster(merge, "predicted.major_cluster")
merge$corrected.minor_cluster <- correct_cluster(merge, "predicted.minor_cluster")

DimPlot(merge, label = T)
DimPlot(merge, reduction = "umap", group.by = "corrected.cluster", label = TRUE, label.size = 3)
DimPlot(merge, reduction = "umap", group.by = "corrected.major_cluster", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(merge, reduction = "umap", group.by = "corrected.minor_cluster", label = TRUE, label.size = 3, repel = TRUE)

saveRDS(merge, file = "result/merged.rds")

## separate stromal, immune and epithelial cells, mapping respectively-----------
map_obj_ref <- function(seurat_obj, reference, resolution = 10){
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  seurat_obj <- seurat_obj %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:40, reduction = "harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  anchors <- FindTransferAnchors(reference = reference, query = seurat_obj, dims = 1:30, 
                                     reference.reduction = "pca")
  
  seurat_obj <- MapQuery(
    anchorset = anchors,
    query = seurat_obj,
    reference = reference,
    refdata = list(
      major_cluster = "major_cluster",
      minor_cluster = "minor_cluster"
    ),
    reference.reduction = "pca", 
    reduction.model = "umap"
  )
  
  seurat_obj$corrected.major_cluster <- correct_cluster(seurat_obj, "predicted.major_cluster")
  seurat_obj$corrected.minor_cluster <- correct_cluster(seurat_obj, "predicted.minor_cluster")
  
  correct_table <- table(seurat_obj$corrected.major_cluster,seurat_obj$corrected.minor_cluster)
  seurat_obj$corrected.major_cluster <- sapply(seurat_obj$corrected.minor_cluster, function(x){
    rownames(correct_table)[which.max(correct_table[,x])]
  })
  
  return(seurat_obj)
}

merged_str <- subset(merge, corrected.cluster == "Stromal")
merged_imm <- subset(merge, corrected.cluster == "Immune")
merged_epi <- subset(merge, corrected.cluster == "Epithelial")

scIBD_subset_str <- subset(scIBD_subset, cluster == "Stromal")
scIBD_subset_imm <- subset(scIBD_subset, cluster == "Immune")
scIBD_subset_epi <- subset(scIBD_subset, cluster == "Epithelial")

merged_str <- map_obj_ref(merged_str, scIBD_subset_str)
merged_imm <- map_obj_ref(merged_imm, scIBD_subset_imm)
merged_epi <- map_obj_ref(merged_epi, scIBD_subset_epi)

DimPlot(merged_str, reduction = "umap", group.by = "corrected.minor_cluster", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(merged_imm, reduction = "umap", group.by = "corrected.minor_cluster", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(merged_epi, reduction = "umap", group.by = "corrected.minor_cluster", label = TRUE, label.size = 3, repel = TRUE)

saveRDS(merged_str, file = "result/merged_str.rds")
saveRDS(merged_imm, file = "result/merged_imm.rds")
saveRDS(merged_epi, file = "result/merged_epi.rds")

## get integrated metadata------
merge_meta <- merged_imm@meta.data[,c("location","group","sample_ident","patient_ident","corrected.cluster","corrected.major_cluster","corrected.minor_cluster")]
scIBD_subset_meta <- scIBD_subset_imm@meta.data[,c("tissue.sub","disease","sample","subject","cluster","major_cluster","minor_cluster")]
colnames(merge_meta) <- colnames(scIBD_subset_meta)

meta_all <- rbind(merge_meta, scIBD_subset_meta)
meta_all$cluster <- as.factor(meta_all$cluster)
meta_all$major_cluster <- as.factor(meta_all$major_cluster)
meta_all$minor_cluster <- as.factor(meta_all$minor_cluster)

meta_all$disease[meta_all$disease=="CEAS"] <- "CEAS_inflamed"
meta_all$disease[meta_all$disease=="SC"] <- "CEAS_non_inflamed"
meta_all$disease[meta_all$disease=="HC"] <- "Healthy"

head(meta_all)
saveRDS(meta_all, file = "result/meta_imm.rds")

## prepare reference data, colon stromal cells-------
stromal_ref <- Read10X(data.dir = "reference/IBD/SCP1884/expression/CO_STR.scp/")
str_ref <- CreateSeuratObject(counts = stromal_ref, project = "stromal_ref", min.cells = 0, min.features = 0)
str_ref

str_meta <- read.table("reference/IBD/SCP1884/metadata/scp_metadata_combined.v2.txt", header = T, row.names = "NAME")
str_meta <- str_meta[-1,]
str_ref <- AddMetaData(object = str_ref, metadata = str_meta)

str_cluster <- read.table("reference/IBD/SCP1884/cluster/CO_STR.scp.X_umap.coords.txt", header = T, row.names = "NAME")
str_cluster <- str_cluster[-1,]
colnames(str_cluster) <- c("UMAP_1","UMAP_2")
str_cluster$UMAP_1 <- as.numeric(str_cluster$UMAP_1)
str_cluster$UMAP_2 <- as.numeric(str_cluster$UMAP_2)
str_cluster <- as.matrix(str_cluster)
umap_reduction <- CreateDimReducObject(embeddings = as.matrix(str_cluster), key = "UMAP_", assay = DefaultAssay(str_ref))
str_ref@reductions$umap <- umap_reduction

predicted_celltype_color_str <- color_database[1:length(unique(str_ref$Celltype))]
names(predicted_celltype_color_str) <- unique(str_ref$Celltype)

pdf("plot_new/basic_CO.str_Crohn_Immunity.pdf")
DimPlot(str_ref, group.by = "Celltype", cols = predicted_celltype_color_str)
DimPlot(str_ref, group.by = "Type")
DimPlot(str_ref, group.by = "disease__ontology_label")
dev.off()

## mapping on reference, colon stromal cells-----------
str_ref <- NormalizeData(str_ref)
str_ref <- FindVariableFeatures(str_ref)
str_ref <- ScaleData(str_ref)
str_ref <- RunPCA(str_ref)

# select two technologies for the query datasets
anchors <- FindTransferAnchors(reference = str_ref, query = merge_str, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = str_ref$Celltype, dims = 1:30)
merge_str <- AddMetaData(merge_str, metadata = predictions)

table(merge_str$predicted.id)
DimPlot(merge_str, group.by = "predicted.id", cols = predicted_celltype_color_str)

## prepare query data, colon and ileum immune cells--------
DimPlot(merge, group.by="celltype")
merge_imm <- subset(merge, celltype %in% c("B_cell","T_cell","NK_cell","Mast_cell", "Neutrophil", "Monocyte/Macrophage") &
                      location %in% c("Colon","Ileum"))
merge_imm <- RunUMAP(merge_imm, dims = 1:30)
DimPlot(merge_imm, group.by = "celltype")

## prepare reference data, , colon and ileum immune cells------
immune_ref_TI <- Read10X(data.dir = "reference/IBD/SCP1884/expression/TI_IMM.scp/")
imm_ref_TI <- CreateSeuratObject(counts = immune_ref, project = "immune_ref", min.cells = 0, min.features = 0)
imm_ref_TI

immune_ref_CO <- Read10X(data.dir = "reference/IBD/SCP1884/expression/CO_IMM.scp/")
imm_ref_CO <- CreateSeuratObject(counts = immune_ref_CO, project = "immune_ref", min.cells = 0, min.features = 0)
imm_ref_CO

imm_ref <- merge(imm_ref_TI, y = list(imm_ref_CO))

imm_meta <- read.table("reference/IBD/SCP1884/metadata/scp_metadata_combined.v2.txt", header = T, row.names = "NAME")
imm_meta <- imm_meta[-1,]
imm_ref <- AddMetaData(object = imm_ref, metadata = imm_meta)

imm_cluster_TI <- read.table("reference/IBD/SCP1884/cluster/TI_IMM.scp.X_umap.coords.txt", header = T, row.names = "NAME")
imm_cluster_TI <- imm_cluster_TI[-1,]

imm_cluster_CO <- read.table("reference/IBD/SCP1884/cluster/CO_IMM.scp.X_umap.coords.txt", header = T, row.names = "NAME")
imm_cluster_CO <- imm_cluster_CO[-1,]

imm_cluster <- rbind(imm_cluster_CO, imm_cluster_TI)
imm_cluster <- imm_cluster[colnames(imm_ref),]

colnames(imm_cluster) <- c("UMAP_1","UMAP_2")
imm_cluster[,1] <- as.numeric(imm_cluster[,1])
imm_cluster[,2] <- as.numeric(imm_cluster[,2])
imm_cluster <- as.matrix(imm_cluster)
umap_reduction <- CreateDimReducObject(embeddings = as.matrix(imm_cluster), key = "UMAP_", assay = DefaultAssay(imm_ref))
imm_ref@reductions$umap <- umap_reduction

predicted_celltype_color_imm <- color_database[(length(unique(imm_ref$Celltype))+1):(length(unique(imm_ref$Celltype))+length(imm_ref$Celltype))]
names(predicted_celltype_color_imm) <- unique(imm_ref$Celltype)

pdf("plot_new/basic_COandTI.imm_Crohn_Immunity.pdf")
DimPlot(imm_ref, group.by = "Celltype", cols = predicted_celltype_color_imm)
DimPlot(imm_ref, group.by = "Type")
DimPlot(imm_ref, group.by = "disease__ontology_label")
dev.off()

## mapping on reference, colon and ileum immune cells-----------
imm_ref <- NormalizeData(imm_ref)
imm_ref <- FindVariableFeatures(imm_ref)
imm_ref <- ScaleData(imm_ref)
imm_ref <- RunPCA(imm_ref)

# select two technologies for the query datasets
anchors <- FindTransferAnchors(reference = imm_ref, query = merge_imm, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = imm_ref$Celltype, dims = 1:30)
merge_imm <- AddMetaData(merge_imm, metadata = predictions)

table(merge_imm$predicted.id)
DimPlot(merge_imm, group.by = "predicted.id", cols = predicted_celltype_color_imm, label = T)

pdf("plot_new/merge_map2ref_imm_and_str.pdf", width = 8, height = 5)
DimPlot(merge_str, group.by = "predicted.id", cols = predicted_celltype_color_str)
DimPlot(merge_imm, group.by = "predicted.id", cols = predicted_celltype_color_imm)
dev.off()

## differential expression of endothelial cells---------
endo_ref <- subset(scIBD_subset, grepl("Endothelial", scIBD_subset$major_cluster) & disease %in% c("CD_inflamed", "UC_inflamed"))
endo_que <- readRDS("result/endothelial_cells_rmbatch.rds")
endo_que <- subset(endo_que, group == "CEAS")
DimPlot(endo_que, group.by = "seurat_clusters")

# set UMAP models
umap_new_model <- list()
umap_new_model$n_epochs <- 500
umap_new_model$alpha <-1
umap_new_model$method <- "umap"
umap_new_model$negative_sample_rate <- 5
umap_new_model$gamma <- 1
umap_new_model$approx_pow <- 0
umap_new_model$n_neighbors <- 30
umap_new_model$metric$cosine <- list()
umap_new_model$embedding <- endo_que[["umap"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap_new_model$a <- ab_param["a"]
umap_new_model$b <- ab_param["b"]
endo_que[["umap"]]@misc$model <- umap_new_model
DimPlot(endo_que, group.by = "seurat_clusters")

anchors <- FindTransferAnchors(reference = endo_que, query = endo_ref, dims = 1:30,
  reference.reduction = "pca")
endo_ref <- MapQuery(
    anchorset = anchors,
    query = endo_ref,
    reference = endo_que,
    refdata = list(
    seurat_clusters = "seurat_clusters"
  ),
    reference.reduction = "pca",
    reduction.model = "umap"
)

DimPlot(endo_que, group.by = "seurat_clusters")
DimPlot(endo_ref, reduction = "ref.umap", group.by = "predicted.seurat_clusters")

# find markers
endo_que$ref_que_label <- "que"
endo_ref$ref_que_label <- "ref"
endo_merge <- merge(endo_que, y = list(endo_ref))

Idents(endo_merge)<-"ref_que_label"
endo_merge <- JoinLayers(endo_merge)
markers <- FindMarkers(endo_merge, ident.1 = "que",ident.2 = "ref", 
                       test.use = "MAST",
                       latent.vars = "ref_que_label")
# Remove ribosomal, mitochondrial, and HLA genes
markers <- markers[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(markers)),]

write.csv(markers, file = "MAST_endo_CEASvsCD_markers.csv", quote = F, col.names = T, row.names = T)
# the volcano plot of endothelial cells-------
library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 5, fdr = 0.0001)
data$row <- rownames(data)
data$regulate[data$regulate=="Up"] <- "CEAS"
data$regulate[data$regulate=="Down"] <- "IBD"
data$regulate[data$regulate=="Normal"] <- "non-significant"

pdf("plot_new/volcano_endo_CEASvsCD.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 5, FDR_cut = 0.0001, legend_position = "DL",
          colors = c("IBD"="navy","CEAS"="firebrick3","non-significant"="grey"), 
          fills = c("IBD"="navy","CEAS"="firebrick3","non-significant"="grey"), 
          label = "row", label_number = 30, output = FALSE)
dev.off()

## GO enrichment of endothelial cells----------
GO_a <- function(interest_gene, global_gene) {  
  require(clusterProfiler)    # Convert gene symbols to ENTREZ IDs for interest and global genes  
  gene_entrez_id2GO <- clusterProfiler::bitr(interest_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID    
  universe_gene_entrez_id2GO <- clusterProfiler::bitr(global_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID    # Perform GO enrichment analysis and simplify results  
  enr_res <- clusterProfiler::enrichGO(    
    gene = gene_entrez_id2GO,    
    universe = universe_gene_entrez_id2GO,    
    ont = "ALL",    
    OrgDb = 'org.Hs.eg.db')  
  enr_res2 <- clusterProfiler::simplify(enr_res)    # Generate plots  
  p1 <- goplot(enr_res)  
  p2 <- barplot(enr_res2, showCategory = 20, color = "p.adjust")    # Return results and plots as a list  
  list(go_res = enr_res,    
       go_res_simplify = enr_res2,    
       p1 = p1,p2 = p2)
  }

GCB_c01 = markers %>% dplyr::filter(p_val_adj < 0.001, avg_log2FC < -3)
GCB_c02 = markers %>% dplyr::filter(p_val_adj < 0.001, avg_log2FC > 3)

CD_go = GO_a(rownames(GCB_c01),rownames(endo_merge))
CEAS_go = GO_a(rownames(GCB_c02),rownames(endo_merge))

process_GO <- function(GO_data, annotation) {  
  GO_data$go_res_simplify@result %>%    
    head(10) %>% 
    mutate(annotation = annotation,
                       log10q = -log10(abs(qvalue)))
}

GO_result <- bind_rows( 
  process_GO(CEAS_go, "CEAS"), 
  process_GO(CD_go, "CD"))

GO_result$log10q = ifelse(GO_result$annotation == "CD", GO_result$log10q * -1, GO_result$log10q)
GO_result$label_hjust=ifelse(GO_result$log10q>0,1,0)

# plot the GO result
pdf("plot_new/GO_endo_CEASvsCD.pdf")
ggplot(GO_result) +  
  geom_bar(aes(    
    x = reorder(Description, log10q),    
    y = log10q,    
    fill = annotation  
    ),  stat = "identity",  color = "white") +  
  geom_text(aes(    
    x = reorder(Description, log10q),    
    y = 0,    
    label = Description,    
    hjust = label_hjust  
    ),  size = 7 * 0.35,  angle = 0) +  
  scale_fill_manual(    
    name = "",    
    values = c("CD" = "#0f7b9f", "CEAS" = "#d83215")  
    ) +  coord_flip() +  cowplot::theme_cowplot() +  
  theme(    
    axis.text.x = element_text(size = 7),    
    text = element_text(size = 7),    
    axis.line = element_line(size = 0.3),    
    axis.ticks = element_line(size = 0.3),    
    axis.text.y = element_blank(),    
    axis.ticks.y = element_blank(),    
    axis.line.y = element_blank(),    
    panel.border = element_blank(),    
    legend.position = "none"  
    ) +  labs(x = "", y = "CD vs CEAS") +  ylim(-10, 12)
dev.off()


