## Figure 5A: volcano plot of endothelial cells----------
endo <- readRDS("result/endothelial_cells_rmbatch.rds")
Idents(endo) <- "group"
markers <- FindMarkers(endo, ident.1 = "CEAS", ident.2 = "HC",
  test.use = "wilcox", logfc.threshold = 0)

markers_filter <- markers %>%
  dplyr::mutate(gene = rownames(markers)) %>%
  filter(!str_starts(gene, "MT-"))

CEAS_down <- markers_filter %>% arrange(avg_log2FC) %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% 
  pull(gene)

pdf("figures/Fig5A.pdf", width = 5, height = 6)
EnhancedVolcano(markers_filter,
  lab = markers_filter$gene,               # Label the genes
  x = 'avg_log2FC',                   # Log fold-change (x-axis)
  y = 'p_val_adj',                    # Adjusted p-value (y-axis)
  title = "Volcano Plot",             # Title of the plot
  pCutoff = 0.01,                     # Significance threshold for p-value
  FCcutoff = 0.6,                    # Fold-change threshold
  pointSize = 1,                      # Size of the points
  labSize = 3,                        # Size of the labels
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.9,
  drawConnectors = TRUE, 
  xlim = c(-3,3),
  axisLabSize = 14,
  hline = c(10e-8)
)
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
GO_result <- enrichGO(gene  = CEAS_down,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable     = TRUE)
dotplot(GO_result)

## Figure 5B: UMAP of endothelial cells-----------
GroupColor <- readRDS("color/GroupColor.rds")

pdf("figures/Fig5B.pdf", width = 6, height = 5)
DimPlot(endo, group.by = "seurat_clusters", label = T, label.size = 3)
DimPlot(endo, group.by = "group", cols = GroupColor)
dev.off()

## Figure 5C: Milo analysis
# Milo is an analysisi for differential abundance
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)

endo <- readRDS("result/endothelial_cells_rmbatch.rds")

traj_sce <- as.SingleCellExperiment(endo)
plotUMAP(traj_sce)

traj_milo <- Milo(traj_sce)
traj_milo <- buildGraph(traj_milo, k = 10, d = 10, reduced.dim="PCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d = 10, reduced_dim="PCA", refinement_scheme = "graph",refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample_ident")

plotNhoodSizeHist(traj_milo)
head(nhoodCounts(traj_milo))

traj_design <- endo@meta.data[,c("sample_ident", "group")]
traj_design <- dplyr::distinct(traj_design) %>% 
  arrange(group) %>% 
  mutate(group = factor(group, levels = c("HC","CEAS")))
rownames(traj_design) <- traj_design$sample_ident

traj_design

da_results <- testNhoods(traj_milo, design = ~ group, design.df = traj_design, fdr.weighting = "graph-overlap")

da_results %>%
  arrange(SpatialFDR) %>%
  head()

traj_milo <- buildNhoodGraph(traj_milo)
umap_pl <- plotReducedDim(traj_milo, dimred = "UMAP", 
                          colour_by="group", 
                          text_size = 3, point_size=0.5) + scale_color_manual(values = GroupColor)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,size_range=c(0.8,3), layout="UMAP",alpha=1)+
  scale_fill_gradient2(low='blue', mid='white', high="red")

pdf("figures/Fig5C-2.pdf", height = 4, width = 9)
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
dev.off()

da_results <- annotateNhoods(traj_milo, 
                             da_results, 
                             coldata_col = "seurat_clusters")
head(da_results)

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha=1)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
pDa


## Figure 5D-H: EC types annotation and monocle2---------
setwd("~/gastroenterology")
endo <- readRDS("~/gastroenterology/result/endothelial_cells_rmbatch.rds")
healthy_ref <- readRDS("~/gastroenterology/reference/endothelial/1_Healthy_Pan-GI_atlas_all_lineages_20241119.rds")

healthy_ref_endo_subset <- subset(healthy_ref, subset = control_vs_disease=="control" & level_1_annot=="Endothelial" & 
                                    (organ_groups=="Large_intestine" | organ_groups=="Small_intestine"))

healthy_ref_endo_subset
saveRDS(healthy_ref_endo_subset, file = "reference/endothelial/healthy_ref_endo_subset.rds")

library(Seurat)
workflow_seurat <- function(seu){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 50)
  seu <- RunUMAP(seu, dims = 1:40)
  seu <- FindNeighbors(seu, dims = 1:40)
  seu <- FindClusters(seu)
  return(seu)
}
healthy_ref_endo_subset <- workflow_seurat(healthy_ref_endo_subset)

DimPlot(healthy_ref_endo_subset, group.by = "level_3_annot",label = T)

healthy_ref_endo_subset_rm_lym <- subset(healthy_ref_endo_subset, subset = level_3_annot!="EC_lymphatic")
endo_rm_lym <- subset(endo, subset = EC_type!="lymphatic")

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
umap_new_model$embedding <- healthy_ref_endo_subset_rm_lym[["umap"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap_new_model$a <- ab_param["a"]
umap_new_model$b <- ab_param["b"]
healthy_ref_endo_subset_rm_lym[["umap"]]@misc$model <- umap_new_model

anchors <- FindTransferAnchors(reference = healthy_ref_endo_subset_rm_lym, query = endo_rm_lym, dims = 1:30, 
                               reference.reduction = "pca")
endo_rm_lym <- MapQuery(
  anchorset = anchors,
  query = endo_rm_lym,
  reference = healthy_ref_endo_subset_rm_lym,
  refdata = list(
    cluster = "cluster",
    level_3_annot = "level_3_annot"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(endo_rm_lym, group.by = "predicted.level_3_annot",label = T)
DimPlot(endo_rm_lym, group.by = "group",label = T)

endo_rm_lym$EC_type <- endo_rm_lym$predicted.level_3_annot
endo_rm_lym$EC_type[endo_rm_lym$seurat_clusters==3] <- "CEAS-specific"
healthy_ref_endo_subset_rm_lym$EC_type <- healthy_ref_endo_subset_rm_lym$level_3_annot

merge_withref <- merge(endo_rm_lym, y = list(healthy_ref_endo_subset_rm_lym))
merge_withref <- workflow_seurat(merge_withref)
merge_withref$group[is.na(merge_withref$group)] <- "REF"

saveRDS(endo_rm_lym, file = "result/endo_rm_lym.rds")

library(harmony)
merge_withref <- merge_withref %>% 
  RunHarmony("group", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
merge_withref <- merge_withref %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 
merge_withref <- merge_withref %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

DimPlot(merge_withref, group.by = "EC_type",label = T)
DimPlot(merge_withref, group.by = "group",label = T)

saveRDS(merge_withref, file = "result/merge_withref_endo.rds")

## markers analysis--------
library(clusterProfiler)

markers <- FindAllMarkers(endo_rm_lym, group.by = "EC_type", only.pos = T)
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
p <- dotplot(xx) + 
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12)) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50)) +
  coord_flip()  # Flip the axes

pdf("figures/Fig5I.pdf", height = 5, width = 10)
print(p)
dev.off()

## pseudotime analysis-------
library(patchwork)
library(monocle)
library(dplyr)

DefaultAssay(merge_withref) <- 'RNA'
cds <- as.CellDataSet(merge_withref)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
## ordering by marker gene per cluster
Idents(merge_withref) <- "EC_type"
deg <- FindAllMarkers(merge_withref)
dim(deg)
deg_filter = deg %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
sel.gene <- unique(deg_filter$gene)

cds <- monocle::setOrderingFilter(cds, sel.gene)

## dimension reduciton
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

## ordering cells
cds <- monocle::orderCells(cds)

library(ggpubr)
library(ggsci)
df<-cds@phenoData@data
GroupColor <- readRDS("color/GroupColor.rds")
cds$group <- factor(cds$group, levels = c("REF","HC","CEAS"))
saveRDS(cds, file="result/merge_with_ref_cds.rds")

pdf("figures/Fig5D-2.pdf", height = 4, width = 5.5)
DimPlot(endo_rm_lym, group.by = "EC_type") + scale_color_nejm()
dev.off()

pdf("figures/Fig5E-F-1.pdf")
monocle::plot_cell_trajectory(cds, color_by = "group", cell_size = 0.5) + scale_color_manual(values= c(GroupColor,"REF"="lightgrey"))
monocle::plot_cell_trajectory(cds, color_by = "EC_type", cell_size = 0.5)+scale_color_nejm()
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1)
dev.off()

pdf("figures/Fig5E-F-2.pdf", height = 3, width = 4.5)
ggplot(df,aes(Pseudotime,color=group,fill=group))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_manual(values=c(GroupColor,"REF"="lightgrey"))+scale_color_manual(values=c(GroupColor,"REF"="lightgrey"))
dev.off()

pdf("figures/Fig5E-F-3.pdf", height = 3, width = 5)
ggplot(df,aes(Pseudotime,color=EC_type,fill=EC_type))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_nejm()+scale_color_nejm()
dev.off()

## Figure 5J: gene module---------
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")

selected_genes <- rownames(subset(diff_test_res,
                    qval < 1e-4))
CEAS_markers <- deg %>% filter(cluster=="CEAS-specific" & p_val_adj < 0.05 & avg_log2FC > 1) %>% pull(gene)
# other_markers <- deg %>% filter(cluster %in% c("EC_capillary","EC_arterial","EC_arterial_2",
#     "EC_cycling","EC_venous") & p_val_adj < 0.01 & avg_log2FC > 1) %>% pull(gene)

filtered_names <- intersect(selected_genes,CEAS_markers)
# filtered_names <- setdiff(selected_genes, other_markers)
length(filtered_names)

cds_subset <- cds[filtered_names, ]

ht <- plot_pseudotime_heatmap(cds_subset,
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = F, 
                        return_heatmap = T)

pdf("figures/Fig5J.pdf", height = 6, width = 5)
ht
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)

gene_modules <- cutree(ht$tree_row, k = 4)
gene_module_df <- data.frame(
  gene_id = names(gene_modules),
  module = as.factor(gene_modules)
)

modules <- split(gene_module_df$gene_id, gene_module_df$module)

library(clusterProfiler)
library(org.Hs.eg.db)
GO_module <- enrichGO(gene  = modules[[1]],
         OrgDb        = org.Hs.eg.db,
         keyType      = "SYMBOL",
         ont          = "BP",
         pAdjustMethod = "BH",
         qvalueCutoff = 0.05,
         readable     = TRUE)
dotplot(GO_module)

## beyondcell
library(beyondcell, lib.loc = "/home/fengqian/miniconda3/envs/beyondcell/lib/R/library")

gs <- GetCollection(PSc)
#ss <- GetCollection(SSc)

sc <- endo_rm_lym

bc <- bcScore(sc, gs, expr.thres = 0.1) 
#bc <- bcScore(ss, ss, expr.thres = 0.1)

# Run the UMAP reduction. 
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.2)

bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 0.5)

bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")
bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

bcRanks <- function(bc, idents = NULL, extended = TRUE, 
                    resm.cutoff = c(0.1, 0.9), 
                    sp.cutoff = c(0.1, 0.4, 0.6, 0.9)) {
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (idents %in% colnames(bc@meta.data)) {
      if (idents %in% names(bc@ranks)) {
        warning(paste0('$', idents, ' already exists in bc@ranks. ',
                       'Entry has been overwritten.'))
      }
      meta <- bc@meta.data[colnames(bc@normalized), idents, drop = TRUE]
    } else {
      stop('Idents not found.')
    }
  } else {
    stop("You must supply the name of a metadata column to group by.")
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # Check resm.cutoff.
  if (length(resm.cutoff) != 2 | !is.numeric(resm.cutoff)) {
    stop('resm.cutoff must be a numeric vector of length 2.')
  }
  if (resm.cutoff[2] < resm.cutoff[1]) {
    warning(paste('Upper residuals\' mean cut-off is smaller than lower', 
                  'residuals\' mean cut-off. Sorting residuals\' mean cut-offs', 
                  'in increasing order.'))
    resm.cutoff <- sort(resm.cutoff, decreasing = FALSE)
  }
  # Check sp.cutoff.
  if (length(sp.cutoff) != 4 | !is.numeric(sp.cutoff)) {
    stop('sp.cutoff must be a numeric vector of length 4.')
  }
  if (any(sp.cutoff < 0 | sp.cutoff > 1)) {
    stop('sp.cutoff must contain 4 switch point values between 0 and 1.')
  }
  sorted.sp.cutoff <- sort(sp.cutoff, decreasing = FALSE)
  if (!identical(sp.cutoff, sorted.sp.cutoff)) {
    warning(paste('Sorting switch point cut-offs in increasing order.'))
    sp.cutoff <- sorted.sp.cutoff
  }
  # --- Code ---
  # Progress bar.
  pb <- txtProgressBar(min = 0, max = 100, style = 3, file = stderr())
  bins <- 10
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # Keep column to group by.
  meta <- bc@meta.data %>%
    tibble::rownames_to_column("cells") %>%
    dplyr::select(cells, all_of(idents)) %>%
    dplyr::rename(group.var := !!idents) %>%
    dplyr::mutate(group.var = factor(group.var)) %>%
    unique()
  lvls <- levels(meta$group.var)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 5)
  # Column to order by.
  order.col <- paste0("rank.", levels(meta$group.var)[1])
  # Final column order.
  if (extended) {
    cols.additional <- c("median", "sd", "variance", "min", "max", "prop.na")
  } else {
    cols.additional <- NULL
  }
  cols.stats <- c("rank", "switch.point", "mean", cols.additional, 
                  "residuals.mean", "group")
  cols.stats.level <- tidyr::expand_grid(lvls, cols.stats) %>%
    dplyr::mutate(col.name = paste(cols.stats, lvls, sep = ".")) %>%
    dplyr::pull(col.name)
  # Get switch points.
  sp <- data.frame(switch.point = bc@switch.point) %>%
    tibble::rownames_to_column("IDs")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 10)
  # Compute long normalized BCS.
  normalized.long <- bc@normalized %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cells") %>%
    tidyr::pivot_longer(cols = all_of(sigs), names_to = "IDs", 
                        values_to = "enrichment", values_drop_na = FALSE)
  # Add grouping information and switch point.
  normalized.long <- normalized.long %>%
    dplyr::inner_join(sp, by = "IDs") %>%
    dplyr::inner_join(meta, by = "cells")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 25)
  # Compute mean BCS and residual's mean per signature.
  stats.long <- normalized.long %>%
    dplyr::group_by(IDs) %>%
    dplyr::mutate(mean = round(mean(enrichment, na.rm = TRUE), digits = 2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(resid = enrichment - mean) %>%
    dplyr::group_by(IDs, group.var) %>%
    dplyr::mutate(residuals.mean = round(mean(resid, na.rm = TRUE), digits = 2)) %>%
    dplyr::ungroup()
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 45)
  # If extended == TRUE, compute the median, standard deviation, variance, min, 
  # max and proportion of NaNs per signature.
  if (extended) {
    stats.long <- stats.long %>%
      dplyr::group_by(IDs) %>%
      dplyr::mutate(median = round(median(enrichment, na.rm = TRUE), digits = 2),
                    sd = round(sd(enrichment, na.rm = TRUE), digits = 2),
                    variance = round(var(enrichment, na.rm = TRUE), digits = 2),
                    min = round(min(enrichment, na.rm = TRUE), digits = 2),
                    max = round(max(enrichment, na.rm = TRUE), digits = 2),
                    prop.na = round(sum(is.na(enrichment))/length(cells), 
                                    digits = 2)) %>%
      dplyr::ungroup()
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 50)
  # Residual's deciles.
  res.decil <- stats.long %>%
    dplyr::group_by(group.var) %>%
    dplyr::group_modify(~as.data.frame(t(quantile(.$residuals.mean, resm.cutoff, 
                                                  na.rm = TRUE)))) %>%
    dplyr::ungroup()
  colnames(res.decil)[2:3] <- c("Pmin", "Pmax")
  stats.long <- stats.long %>%
    dplyr::select(-cells, -enrichment) %>%
    unique() %>%
    dplyr::inner_join(res.decil, by = "group.var")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 75)
  # Group annotation.
  stats.long.annotated <- stats.long %>%
    dplyr::mutate(
      group = dplyr::case_when(switch.point < sp.cutoff[1] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-HighSensitivity",
                               switch.point > sp.cutoff[4] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-LowSensitivity",
                               switch.point > sp.cutoff[2] & 
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-Differential-LowSensitivity",
                               switch.point > sp.cutoff[2] &
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-Differential-HighSensitivity",
                               TRUE ~ NA_character_))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 80)
  # Order.
  rank <- stats.long.annotated %>%
    dplyr::mutate(in.range = switch.point > sp.cutoff[2] & 
                    switch.point < sp.cutoff[3],
                  sp.rank = switch.point * as.numeric(in.range)) %>%
    dplyr::select(IDs, group.var, sp.rank, residuals.mean, in.range) %>%
    unique() %>%
    dplyr::group_split(group.var)
  rank <- lapply(rank, FUN = function(x) {
    dt <- data.table::as.data.table(x)
    dt[, rank := data.table::frank(dt, -sp.rank, -residuals.mean, 
                                   ties.method = "dense")]
    return(dt)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(rank = dplyr::if_else(in.range, rank, NA_integer_)) %>%
    dplyr::select(IDs, group.var, rank) %>%
    unique()
  stats.long.ranked <- stats.long.annotated %>%
    dplyr::inner_join(rank, by = c("IDs", "group.var"))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 85)
  # Pivot wider
  final.stats <- stats.long.ranked %>%
    dplyr::select(IDs, group.var, all_of(cols.stats)) %>%
    unique() %>%
    tidyr::pivot_wider(names_from = group.var, values_from = all_of(cols.stats),
                       names_sep = ".")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 90)
  # Add Drug name and MoA to final.stats.
  info <- drugInfo$IDs %>%
    dplyr::filter(IDs %in% final.stats$IDs) %>%
    dplyr::select(IDs, preferred.drug.names, studies) %>%
    dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Targets, by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Synonyms, by = "IDs",
                     relationship = "many-to-many")
  if (dim(info)[1] > 0) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(x) {
      paste(na.omit(unique(x)), collapse = "; ")
    })
    cols.druginfo <- c("drugs", "preferred.drug.names", "MoAs", "targets", 
                       "studies")
  } else {
    info <- data.frame(IDs = rownames(bc@normalized))
    cols.druginfo <- NULL
  }
  final.stats <- final.stats %>%
    dplyr::left_join(info, by = "IDs") %>%
    tibble::column_to_rownames("IDs") %>%
    unique()
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 95)
  # Order by rank and reorder columns.
  final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                             c(cols.druginfo, cols.stats.level)]
  # Add to beyondcell object.
  bc@ranks[[idents]] <- final.stats
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 100)
  return(bc)
}

bc <- bcRanks(bc, idents = "EC_type")

# Step 1: Extract data from bc@ranks
df <- bc@ranks$EC_type

df <- df %>%
  dplyr::select(
    residuals.mean = "residuals.mean.CEAS-specific",
    switch.point = "switch.point.CEAS-specific",
    group = "group.CEAS-specific",
    drug = preferred.drug.names,
    description = MoAs
  ) %>%
  tibble::rownames_to_column("sig_id") %>%
  dplyr::mutate(
    label = paste0(drug, " (", sig_id, ")")
  )

# Step 2: Label only top 5 from each group
df <- df %>%
  arrange(group, desc(abs(residuals.mean))) %>%
  group_by(group) %>%
  mutate(label_show = ifelse(
    (group == "TOP-Differential-HighSensitivity" & row_number() <= 10) | 
      (group == "TOP-Differential-LowSensitivity" & row_number() <= 3),
    label, 
    NA
  )) %>%
  ungroup()

# Step 3: Define color map
color_map <- c(
  "TOP-Differential-HighSensitivity" = "#F4A300",
  "TOP-Differential-LowSensitivity" = "#B086CC"
)

# Step 4: Set cutoffs
x_cutoff <- quantile(df$residuals.mean, probs = c(0.1, 0.9), na.rm = TRUE)
y_cutoff <- c(0.1, 0.4, 0.6, 0.9)

# Step 5: Plot
pdf("figures/Fig5K.pdf", height = 5, width = 8)
ggplot(df, aes(x = residuals.mean, y = switch.point)) +
  geom_point(aes(fill = group, color = group), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color_map, na.value = "lightgrey") +
  scale_color_manual(values = color_map, na.value = "lightgrey", guide = "none") +
  geom_vline(xintercept = x_cutoff, linetype = "dotted") +
  geom_hline(yintercept = y_cutoff, linetype = "dotted") +
  geom_text_repel(
    aes(label = label_show),
    size = 3,
    color = "black",
    nudge_x = ifelse(df$residuals.mean >= 0, 0.5, -0.5),
    direction = "y",      
    hjust = 0.5,           
    segment.size = 0.5,   
    segment.alpha = 0.6,  
    min.segment.length = 0, 
    box.padding = 0.5,  
    point.padding = 0.3,  
    force = 2,       
    max.time = 2,      
    max.iter = 20000,  
    max.overlaps = Inf   
  ) +
  theme_classic() +
  labs(
    title = "EC_type = CEAS-specific",
    x = "Residuals' Mean",
    y = "Switch Point",
    fill = "Drug Group",
    caption = paste0("x cut-offs: first and last deciles; y cut-offs: ", paste(y_cutoff, collapse = ", "))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
dev.off()

write.csv(df, file = "tables/Sup_Table2.csv", quote = F)

pdf("figures/Sup_Fig5E.pdf")
bcHistogram(bc, signatures = "sig-1439", idents = "EC_type")
bcHistogram(bc, signatures = "sig-20560", idents = "EC_type")
bcHistogram(bc, signatures = "sig-9639", idents = "EC_type")
dev.off()







