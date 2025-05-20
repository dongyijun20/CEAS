# Figure 4A/B: stromal and epithelial cells umap with minor cluster-------
merged_str = readRDS("result/merged_str_lower_GI_tract.rds")
merged_epi = readRDS("result/merged_epi_lower_GI_tract.rds")

color_database = readRDS("color/color_database.rds")
color_str <- color_database[22:(22+length(unique(merged_str$corrected.minor_cluster)))]
names(color_str) <- unique(merged_str$corrected.minor_cluster)

p1 <- DimPlot(merged_str, reduction = "umap", group.by = "corrected.minor_cluster", 
        label = TRUE, label.size = 3, repel = TRUE, cols = color_str) +
  ggtitle("CEAS Stromal Cells UMAP")
ggsave("figures/Fig4A.svg", plot = p1, width = 7, height = 5, dpi = 300)


color_database2 = c('#d4de9c','#94c58f','#9cd2ed','#a992c0',
                    '#ea9994','#f2c396','#bb82b1')
color_epi <- color_database2[1:length(unique(merged_epi$corrected.minor_cluster))]
names(color_epi) <- unique(merged_epi$corrected.minor_cluster)

p2 <- DimPlot(merged_epi, reduction = "umap", group.by = "corrected.minor_cluster", 
        label = TRUE, label.size = 3, repel = TRUE, cols = color_epi) +
  ggtitle("CEAS Epithelial Cells UMAP") 
ggsave("figures/Fig4B.svg", plot = p2, width = 6.5, height = 5, dpi = 300)

## Figure 4C: dot plot for epithelial and mesenchymal -------------
library(RColorBrewer)

prepare_markers_gene <- function(Epi_subset, celltype, top = 5){
  de_results <- list()
  clusters <- unique(Epi_subset$corrected.minor_cluster)
  
  for (cluster in clusters) {
    cluster_cells <- subset(Epi_subset, corrected.minor_cluster == cluster)
    
    de_results[[cluster]] <- FindMarkers(
      cluster_cells, 
      ident.1 = "Inflamed", 
      ident.2 = "Healthy", 
      group.by = "disease", 
      min.pct = 0.1, 
      logfc.threshold = 0.25
    )
  }
  saveRDS(de_results, file = paste0("result/CEASvsHC_markers_",celltype,".rds"))
  
  Idents(Epi_subset) <- "corrected.minor_cluster"
  all_markers <- FindAllMarkers(
    object = Epi_subset,
    only.pos = TRUE,            # Only positive markers
    min.pct = 0.1,              # Minimum expression percentage threshold
    logfc.threshold = 0.25      # Log fold-change threshold
  )
  saveRDS(all_markers, file = paste0("result/all_markers_",celltype,".rds"))
  
  top_genes <- lapply(names(de_results), function(name) {
    if(nrow(de_results[[name]])>0){
      df <- de_results[[name]] %>% 
        arrange(p_val_adj, abs(avg_log2FC)) %>%  
        mutate(gene = rownames(de_results[[name]])) %>% 
        filter(gene %in% all_markers$gene[all_markers$cluster==name & all_markers$p_val_adj < 0.01]) %>% 
        filter(!str_starts(gene, "MT-") & !str_starts(gene, "LINC") & !str_detect(gene, ".*\\.[0-9]+$")) %>% 
        head(top) %>%             # Select top 10 genes (adjust as needed)
        rownames()
      return(df)
    } else return(NULL)
  })
  
  # Flatten the list to a unique vector of gene names
  top_genes <- unique(unlist(top_genes))
  
  top_genes
  
  # Create a new column combining minor cluster and disease condition
  Epi_subset$combined_cluster <- factor(paste(Epi_subset$corrected.minor_cluster, Epi_subset$disease, sep = "_"),
                                        levels = paste(rep(names(de_results),each = 2), c("Healthy","Inflamed"), sep = "_"))
  Idents(Epi_subset) <- "combined_cluster"
  
  return(list(Epi_subset, top_genes))
}

Epi_subset <- subset(merged_epi, disease %in% c("Inflamed","Healthy"))
Mes_subset <- subset(merged_str, corrected.major_cluster == "Mesenchymal" & disease %in% c("Inflamed","Healthy"))

Epi_subset_list <- prepare_markers_gene(Epi_subset, "Epithelial")
Mes_subset_list <- prepare_markers_gene(Mes_subset, "Mesenchymal", top = 6)

# Create DotPlot for the selected top DEGs across clusters and disease conditions
pdf("figures/Fig4E-1.pdf", height = 3.5, width = 9)
DotPlot(
  Epi_subset_list[[1]], 
  features = Epi_subset_list[[2]], 
  group.by = "combined_cluster",  # Use combined cluster (cluster + disease) for separation
  cols = c("white", "coral"),     # Color gradient for expression levels
  dot.scale = 6                   # Adjust dot size
) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradient(low = "grey", high = "coral") +        # Color scale for mean expression
  labs(y = NULL, x = NULL)
dev.off()

pdf("figures/Fig4E-2.pdf", height = 2.5, width = 10)
DotPlot(
  Mes_subset_list[[1]], 
  features = Mes_subset_list[[2]], 
  group.by = "combined_cluster",  # Use combined cluster (cluster + disease) for separation
  cols = c("white", "purple"),     # Color gradient for expression levels
  dot.scale = 6 ,                   # Adjust dot size
) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradient(low = "grey", high = "purple") +        # Color scale for mean expression
  labs(y = NULL, x = NULL)
dev.off()

## function analysis
library(clusterProfiler)

de_results <- readRDS("result/CEASvsHC_markers_Epithelial.rds")
sce.markers <- lapply(names(de_results), function(x){
  de_results[[x]] %>% 
    mutate(celltype = x, gene = rownames(de_results[[x]])) %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0)
})
sce.markers <- Reduce(bind_rows, sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

gcSample=split(sce.markers$ENTREZID, sce.markers$celltype) 
compare1 <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01
)

plot1 <- dotplot(compare1, showCategory=5) + 
  scale_fill_gradientn(name = "p.adjust", colours = brewer.pal(n = 9, name = "RdYlBu")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  ggtitle("CEAS vs HC upregulated genes in Epithelial")
plot1$layers[[1]]$aes_params$stroke <- 0

pdf("figures/Fig4F-1.pdf", height = 6.5, width = 5.5)
plot1
dev.off()

de_results <- readRDS("result/CEASvsHC_markers_Mesenchymal.rds")
sce.markers <- lapply(names(de_results), function(x){
  de_results[[x]] %>% 
    mutate(celltype = x, gene = rownames(de_results[[x]])) %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0)
})
sce.markers <- Reduce(bind_rows, sce.markers)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

gcSample=split(sce.markers$ENTREZID, sce.markers$celltype) 
compare2 <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01
)

plot2 <- dotplot(compare2, showCategory=5) + 
  scale_fill_gradientn(name = "p.adjust", colours = brewer.pal(n = 9, name = "RdYlBu")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("CEAS vs HC upregulated genes in Mesenchymal")
plot2$layers[[1]]$aes_params$stroke <- 0

pdf("figures/Fig4F-2.pdf", height = 7, width = 5.5)
plot2
dev.off()

## Figure4G: Nichenet Cell-cell commmunication---------
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)

## CEAS, endothelial cells are receiver--------
merged_str <- readRDS("result/merged_str_lower_GI_tract.rds")
merged_imm <- readRDS("result/merged_imm_lower_GI_tract.rds")
merged_epi <- readRDS("result/merged_epi_lower_GI_tract.rds")
merge <- readRDS("result/merged_lower_GI_tract.rds")

change_group <- function(merged_str){
  merged_str$disease <- merged_str$group
  merged_str$disease[merged_str$disease=="CEAS"] <- "Inflamed"
  merged_str$disease[merged_str$disease=="SC"] <- "Non_inflamed"
  merged_str$disease[merged_str$disease=="HC"] <- "Healthy"
  return(merged_str)
}
merged_str <- change_group(merged_str)
merged_imm <- change_group(merged_imm)
merged_epi <- change_group(merged_epi)

seuratObj <- merge(x=merged_str, y=list(merged_imm, merged_epi))
seuratObj@reductions <- merge@reductions

Idents(seuratObj) <- "corrected.major_cluster"

seuratObj@meta.data %>% head()

seuratObj <- alias_to_symbol_seurat(seuratObj, "human")

# Note that the number of cells of some cell types is very low and should preferably be higher for a real application
seuratObj@meta.data$corrected.major_cluster %>% table() 

# B_Plasma        CD4T        CD8T Endothelial  Epithelial         ILC Mesenchymal     Myeloid      Neural 
# 19723        8551        9678        2178       15036         545        6360        6519         661 

DimPlot(seuratObj, reduction = "umap")

organism <- "human"

lr_network <- readRDS("reference/nichenet/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("reference/nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("reference/nichenet/weighted_networks_nsga2r_final.rds")

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] 

head(weighted_networks$lr_sig) 

head(weighted_networks$gr)

receiver = "Endothelial"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("B_Plasma", "CD4T", "CD8T", "Epithelial", "ILC", "Mesenchymal", "Myeloid", "Neural")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

condition_oi <-  "Inflamed"
condition_reference <- "Healthy"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "disease",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% dplyr::arrange(-aupr_corrected) %>% mutate(rank = rank(-aupr_corrected))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% as.data.frame() %>% filter(test_ligand %in% best_upstream_ligands) %>%
  arrange(aupr_corrected) %>% column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 40) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

nrow(active_ligand_target_links)
head(active_ligand_target_links)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 40) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

vis_ligand_receptor_network <- t(vis_ligand_receptor_network[,order_ligands])

p_ligand_receptor <- make_heatmap_ggplot(vis_ligand_receptor_network, 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

best_upstream_ligands_all %in% rownames(seuratObj) %>% table()


# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seuratObj, corrected.major_cluster %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot


celltype_order <- levels(Idents(seuratObj)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "disease",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "corrected.major_cluster",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%
  purrr::reduce(full_join, by = "gene") %>%
  column_to_rownames("gene")

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")

pdf("figures/Fig4G.pdf", width = 20, height = 10)
combined_plot
dev.off()


