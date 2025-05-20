library(mergeect, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(ggplot2)
library(patchwork)
library(ArchR)

## Fig2A/B. plot CEAS data basic information---------
merge <- readRDS("result/merge_location_in_lower_GI_tract.rds")

location_color <- readRDS("color/LocationColor.rds")
Idents(merge) <- "seurat_clusters"

celltype_color<-ArchRPalettes$stallion[1:length(unique(merge$corrected.major_cluster))]
names(celltype_color) <- unique(merge$corrected.major_cluster)

# Generate individual DimPlot objects
p1 <- DimPlot(merge, label = TRUE, label.size = 3) + NoLegend() + ggtitle("")
p2 <- DimPlot(merge, reduction = "umap", group.by = "corrected.cluster", label = TRUE, label.size = 4) + NoLegend() + ggtitle("")
p3 <- DimPlot(merge, reduction = "umap", group.by = "location", cols = location_color) + ggtitle("")
p4 <- DimPlot(merge, reduction = "umap", group.by = "corrected.major_cluster", label = TRUE, label.size = 4, repel = TRUE, cols = celltype_color) + NoLegend() + ggtitle("")

ggsave("figures/Fig2A-1.png", plot = p2, width = 5.5, height = 5, dpi = 300)
ggsave("figures/Fig2A-2.png", plot = p3, width = 5.5, height = 5, dpi = 300)
ggsave("figures/Fig2A-3.png", plot = p4, width = 5.5, height = 5, dpi = 300)

## Fig2C. heatmap regarding markers of each celltype--------
merge_subset <- merge[,sample(colnames(merge),2000)]

Idents(merge_subset) <- "corrected.major_cluster"
markers <- FindAllMarkers(merge_subset, only.pos = T, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
markers_df <- markers_df[-grep("^MT",markers_df$gene),]
write.table(markers_df,file="result/markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

top20.markers = markers_df %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
set.seed(12)
CellInfoS <- merge_subset@meta.data
merge.sce=as.SingleCellExperiment(merge_subset)
library(SingleCellExperiment)
plot.data<-as.data.frame(assay(merge.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)

column_annot <-CellInfoS[,c("corrected.major_cluster","disease"),drop=F]
column_annot$disease = as.factor(as.character(column_annot$disease))
column_annot=with(column_annot, column_annot[order(disease), , drop=F])
column_annot=with(column_annot, column_annot[order(corrected.major_cluster), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]

SampleColors<-readRDS("color/GroupColor.rds")
names(SampleColors) <- c("Inflamed","Healthy","Non_inflamed")
AssignmentColors<-readRDS("color/celltype_color.rds")

column.colors=list()
column.colors[["disease"]]<-SampleColors
column.colors[["corrected.major_cluster"]]<-AssignmentColors
sample=as.matrix(column_annot[,c("disease"),drop=F])
seurat_clusters=as.matrix(column_annot[,c("corrected.major_cluster"),drop=F])

library(ComplexHeatmap)
colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
genes= top2$gene
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows, row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$disease)),title="disease",legend_gp = gpar(fill=SampleColors[levels(as.factor(column_annot$disease))],fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$corrected.major_cluster)),title="celltype",legend_gp = gpar(fill=AssignmentColors[levels(as.factor(column_annot$corrected.major_cluster))],fontsize=5))

pdf("figures/Fig2C.pdf",width=7.5,height=9)
draw(HM,heatmap_legend_list = list(lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()

## Fig2D. proportion difference-------------
library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- merge@meta.data
P1=CellInfo %>% ggplot(aes(x=disease, fill=fct_rev(corrected.major_cluster))) +scale_fill_manual(values = celltype_color)+
  geom_bar(color="black",position = "fill",width = 0.7,size = 0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13),legend.title = element_text(color="black",size=13),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(corrected.major_cluster , fill=corrected.major_cluster))+
  geom_bar(stat="count",colour = "black",width = 0.7,size = 0.3)+  
  scale_fill_manual(values = celltype_color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5), axis.title=element_text(size=15),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)
pdf("figures/Fig2D.pdf",width=8,height=4)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

## Fig2E. SLCO2A1 expression difference-----------
library(rstatix)
library(ggpubr)

merge_subset <- subset(merge, disease %in% c("Inflamed","Healthy"))

plot_df <- data.frame(expression=merge_subset@assays$RNA@data["SLCO2A1",],
                      disease=factor(merge_subset$disease),
                      celltype=factor(merge_subset$corrected.major_cluster))
stat.test <- plot_df %>%
  group_by(celltype) %>%
  wilcox_test(expression ~ disease) %>% # use wilcox
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp<-ggboxplot(
  plot_df, x = "celltype", y = "expression", 
  color = "disease", palette = SampleColors) +
  ggtitle("SLCO2A1 Expression in All Cells") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf("figures/Fig2E.pdf", height = 6, width = 8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0
)
dev.off()

## Fig2F. SLCO2A1 expressed cell proportion difference-----------
merge_subset$SLCO2A1_exp <- ifelse(test = merge_subset@assays$RNA@counts["SLCO2A1",]>0,yes = "1",no = "0")
table(merge_subset$SLCO2A1_exp)
SLCO2A1 <- subset(merge_subset, SLCO2A1_exp=="1")

pro_CEAS<-table(SLCO2A1[,SLCO2A1$disease=="Inflamed"]$corrected.major_cluster)/table(merge_subset[,merge_subset$disease=="Inflamed"]$corrected.major_cluster)
a<-table(SLCO2A1[,SLCO2A1$disease=="Healthy"]$corrected.major_cluster)
a<-a[names(pro_CEAS)]
pro_HC<-a/table(merge_subset[,merge_subset$disease=="Healthy"]$corrected.major_cluster)
df <- data.frame(value=c(pro_CEAS,pro_HC),disease=c(rep("Inflamed",9),rep("Healthy",9)),celltype=c(names(pro_CEAS), names(pro_HC)))

pdf("figures/Fig2F.pdf",width=8,height=4)
ggplot(df, aes(celltype, value, fill=disease))+geom_col(position="dodge")+  
  scale_fill_manual(values = SampleColors)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

## Fig2G. KEGG of three large celltypes, compare difference between inflamed and healthy cells---------------------
# Load required packages
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

## compare all major cell type, up and down regulated pathways
merged_imm = readRDS("result/merged_imm_lower_GI_tract.rds")
merged_str = readRDS("result/merged_str_lower_GI_tract.rds")
merged_epi = readRDS("result/merged_epi_lower_GI_tract.rds")

Idents(merged_imm) <- "disease"
Idents(merged_str) <- "disease"
Idents(merged_epi) <- "disease"
markers_imm <- FindMarkers(merged_imm, ident.1 = "Inflamed",ident.2 = "Healthy", test.use = "wilcox")
markers_str <- FindMarkers(merged_str, ident.1 = "Inflamed",ident.2 = "Healthy", test.use = "wilcox")
markers_epi <- FindMarkers(merged_epi, ident.1 = "Inflamed",ident.2 = "Healthy", test.use = "wilcox")

markers_imm <- markers_imm %>%
  dplyr::mutate(gene = rownames(markers_imm),
                geneID = bitr(rownames(markers_imm), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = F)[,2])
markers_str <- markers_str %>%
  dplyr::mutate(gene = rownames(markers_str),
                geneID = bitr(rownames(markers_str), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = F)[,2])

bitr_results <- bitr(rownames(markers_epi), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
bitr_results <- bitr_results[!duplicated(bitr_results$SYMBOL), ] # Keep only the first mapping

markers_epi <- markers_epi %>%
  dplyr::mutate(
    gene = rownames(markers_epi),
    geneID = bitr_results$ENTREZID[match(rownames(markers_epi), bitr_results$SYMBOL)]
  )

# Function to perform KEGG enrichment analysis
run_kegg <- function(marker_df, cell_type) {
  # Filter significant genes (adjust p-value < 0.05)
  pos_genes <- marker_df %>% 
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    pull(geneID)
  neg_genes <- marker_df %>% 
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
    pull(geneID)
  
  # Run KEGG enrichment analysis
  kegg_results_pos <- as.data.frame(enrichKEGG(gene = na.omit(pos_genes),
                                               organism = 'hsa', 
                                               pvalueCutoff = 0.05)) %>%
    mutate(cell_type = cell_type,
           regulation = "positive")
  kegg_results_neg <- as.data.frame(enrichKEGG(gene = na.omit(neg_genes),
                                               organism = 'hsa', 
                                               pvalueCutoff = 0.05)) %>%
    mutate(cell_type = cell_type,
           regulation = "negative",
           Count = -Count)
  
  # Add cell type information
  kegg_results_df <- bind_rows(kegg_results_pos, kegg_results_neg) %>% 
    arrange(p.adjust)
  
  return(kegg_results_df)
}

# Run KEGG analysis for each cell type
immune_kegg <- run_kegg(markers_imm, "Immune")
epithelial_kegg <- run_kegg(markers_epi, "Epithelial")
stromal_kegg <- run_kegg(markers_str, "Stromal")

# select pathways manually (exclude the unrelated diseases), retain 5 positive and 5 negative terms
# Define keywords to exclude
exclude_keywords <- c("COVID", "Leishmaniasis", "Chagas", "Measles", "Tuberculosis", "Leishmaniasis", "diabet", "thyroid", "Axon", "depression", "Influenza",
                      "Toxoplasmosis", "virus", "addiction", "Hematopoietic", "Graft", "myocarditis", "Allograft", "cancer", "atherosclerosis", "Bile", "cardi",
                      "noma", "leukemia", "Hepatitis", "Yersinia", "Shigellosis", "coli", "Malaria", "Cardiac", "Osteoclast")

# Helper function to filter KEGG results
filter_kegg <- function(kegg_df) {
  kegg_df %>%
    filter(!str_detect(Description, paste(exclude_keywords, collapse = "|"))) %>%
    filter(qvalue < 0.01) %>% 
    arrange(Count) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
}

# Apply filtering to each cell type's KEGG results
immune_kegg_filtered <- filter_kegg(immune_kegg)
epithelial_kegg_filtered <- filter_kegg(epithelial_kegg)
stromal_kegg_filtered <- filter_kegg(stromal_kegg)

plot_kegg <- function(kegg_filtered, cell_type){
  ggplot(kegg_filtered, aes(
    x = Count, 
    y = Description, 
    fill = regulation
  )) +
    geom_col() +
    scale_fill_manual(values = c("positive" = "coral", "negative" = "skyblue")) +
    theme_minimal() +
    labs(
      x = "Counts",
      y = "KEGG Pathways",
      fill = "Regulation",
      title = paste0("KEGG Pathway Enrichment for ", cell_type, " (q-value < 0.01)")
    ) +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

pdf("figures/Fig2G-1.pdf", width = 6, height = 4)
plot_kegg(immune_kegg_filtered, "Immune")
dev.off()

pdf("figures/Fig2G-2.pdf", width = 5, height = 4)
plot_kegg(epithelial_kegg_filtered, "Epithelial")
dev.off()

pdf("figures/Fig2G-3.pdf", width = 6, height = 4)
plot_kegg(stromal_kegg_filtered, "Stromal")
dev.off()

