# Figure 3A: immune cells umap with minor cluster-------
merged_imm = readRDS("result/merged_imm_lower_GI_tract.rds")

color_database = readRDS("color/color_database.rds")
color_imm <- color_database[1:length(unique(merged_imm$corrected.minor_cluster))]
names(color_imm) <- unique(merged_imm$corrected.minor_cluster)

p1 <- DimPlot(merged_imm, reduction = "umap", group.by = "corrected.minor_cluster", 
        label = TRUE, label.size = 3, repel = TRUE, cols = color_imm) +
  ggtitle("CEAS Immune Cells UMAP")

ggsave("figures/Fig3A.png", plot = p1, width = 7, height = 5, dpi = 300)

# Figure 3B: compare different disease-------
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(circlize)

# load meta data
meta = readRDS("result/meta_imm_lower_GI_tract.rds")

# order subset
subset_order = c(meta[meta$major_cluster == "CD4T",]$minor_cluster %>% droplevels() %>% levels,
                 meta[meta$major_cluster == "CD8T",]$minor_cluster %>% droplevels() %>% levels,
                 meta[meta$major_cluster == "ILC",]$minor_cluster %>% droplevels() %>% levels,
                 meta[meta$major_cluster == "B_Plasma",]$minor_cluster %>% droplevels() %>% levels,
                 meta[meta$major_cluster == "Myeloid",]$minor_cluster %>% droplevels() %>% levels)

data.sub = meta %>%
  dplyr::filter(major_cluster %in% c("CD4T","CD8T","ILC","B_Plasma","Myeloid")) %>% 
  dplyr::filter(!grepl("non_inflamed",disease))
  
# count
N.o <- table(data.sub[["minor_cluster"]],data.sub[["disease"]])

# order disease
N.o = N.o[,c("CEAS_inflamed","CD_inflamed","UC_inflamed","Healthy")]
colnames(N.o) <- c("CEAS", "CD", "UC", "HC")

# filter cluster
N.o = N.o[which(rowSums(N.o > 20) > 0) %>% names(),]

# order cluster
N.o = N.o[subset_order[subset_order %in% rownames(N.o)],]

# calculate expected value
res.chisq <- chisq.test(N.o)

# calculate Ro/e
R.oe <- (res.chisq$observed)/(res.chisq$expected)

mark = R.oe
mark[R.oe == 0] = '-'
mark[R.oe > 0 & R.oe < 0.5] = '+/-'
mark[R.oe >= 0.5 & R.oe <= 1] = '+'
mark[R.oe > 1 & R.oe <= 2] = '++'
mark[R.oe > 2] = '+++'

cell_function = function(j, i, x, y, width, height, fill){
  grid.text(mark[i, j], x, y, gp = gpar(fontsize = 8))
}

#getPalette = colorRampPalette(brewer.pal(11, "Spectral")[2:6] %>% rev)
getPalette = colorRampPalette(brewer.pal(9, "Blues")[1:6])
R.oe.v2 = R.oe
R.oe.v2[R.oe.v2 > 3] = 3

row_split = factor(c(rep("1-CD4+ T", unique(meta$minor_cluster[meta$major_cluster=="CD4T" & meta$minor_cluster %in% rownames(R.oe.v2)]) %>% length), 
              rep("2-CD8+ T", unique(meta$minor_cluster[meta$major_cluster=="CD8T" & meta$minor_cluster %in% rownames(R.oe.v2)])  %>% length),
              rep("3-ILC", unique(meta$minor_cluster[meta$major_cluster=="ILC" & meta$minor_cluster %in% rownames(R.oe.v2)])  %>% length), 
              rep("4-B/Plasma",unique(meta$minor_cluster[meta$major_cluster=="B_Plasma" & meta$minor_cluster %in% rownames(R.oe.v2)])  %>% length),
              rep("5-Myeloid", unique(meta$minor_cluster[meta$major_cluster=="Myeloid" & meta$minor_cluster %in% rownames(R.oe.v2)]) %>% length )),
  levels = c("1-CD4+ T", "2-CD8+ T", "3-ILC", "4-B/Plasma", "5-Myeloid"))
  
pdf("figures/Fig3B.pdf", width = 4, height = 10)
Heatmap(R.oe.v2,
        column_title = "",
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_order = factor(c("CEAS", "CD", "UC", "HC"), levels = c("CEAS", "CD", "UC", "HC")),
        row_split = row_split,
        cell_fun = cell_function,
        col = getPalette(5), 
        row_names_gp = gpar(fontsize = 10), 
        border = TRUE,
        heatmap_legend_param = list(
          title = "Ro/e",
          break_dist = 1,
          #col_fun = colorRamp2(c(0, 1, 2, 4, max(R.oe.v2)), getPalette(5)),
          #at = c(0, 1, 2, 4, max(R.oe.v2)),
          col_fun = colorRamp2(c(0, 0.5, 1, 2, 3), getPalette(5)),
          at = c(0, 1, 2, 3, 4),
          labels = c('0','0.5','1','2','max')
          #labels = c('0','0.5','1','4',sprintf("%.1f", max(R.oe.v2)))
        )
)
dev.off()

## Figure 3C: boxplot show different percentage in each sample-------
library(patchwork)

df <- meta %>%
  dplyr::filter(major_cluster %in% c("CD4T","CD8T","ILC","B_Plasma","Myeloid")) %>% 
  dplyr::filter(!grepl("non_inflamed",disease)) %>%
  group_by(subject, major_cluster, minor_cluster, disease) %>%  # Group by subject, major_cluster, and minor_cluster
  add_count(name = "count") %>%     
  dplyr::distinct() %>% 
  group_by(subject, major_cluster, disease) %>%                  # Group by subject and major_cluster to calculate total count
  mutate(total_count = sum(count)) %>%                  # Calculate the total count of cells in each major_cluster
  ungroup() %>%
  mutate(percentage = (count / total_count) * 100) %>%     
  dplyr::distinct() %>%  # Calculate percentage of each minor_cluster within the major_cluster
  select(subject, minor_cluster, disease, percentage)          # Select the relevant columns

df_filtered <- df %>%
  group_by(minor_cluster, disease) %>%
  filter(minor_cluster %in% unique(merged_imm$corrected.minor_cluster) & minor_cluster != "Cycling GC B") %>%
  ungroup() %>% 
  dplyr::distinct()

palette <- c(
  "CD" = "#4DB3E6", 
  "HC" = "#E69F00", 
  "UC" = "#00BA38", 
  "CEAS" = "#F0E442"
)

plot_list <- list()
for(cell_type in unique(df_filtered$minor_cluster)){
  df_filtered_cell_type <- df_filtered %>% 
    filter(minor_cluster == cell_type) %>%
    mutate(disease = recode(disease, 
                            "CEAS_inflamed" = "CEAS",
                            "Healthy" = "HC",
                            "UC_inflamed" = "UC",
                            "CD_inflamed" = "CD"))
  df_filtered_cell_type$disease <- factor(df_filtered_cell_type$disease, levels = c("CEAS","CD","UC","HC"))
  
  stat.test <- df_filtered_cell_type %>%
    wilcox_test(percentage ~ disease) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test <- stat.test %>%
    add_xy_position(x = "disease", dodge = 0.8)
  
  if(sum(stat.test$p.adj.signif!="ns")>0){
    bxp <- ggboxplot(
      df_filtered_cell_type, 
      x = "disease", 
      y = "percentage", 
      fill = "disease", color = "disease",
      palette = palette, alpha = 0.7,
      outlier.shape = NA # Hide outliers to show all points within the box
    ) +
      geom_point(size = 2, aes(col = disease)) +
      scale_color_manual(values = palette) +
      ggtitle(cell_type) +
      theme(
        legend.position = "none" # Remove legend if not needed
      ) +
      labs(y = "Percentage of cells", x = NULL)
    
    plot_list[[cell_type]] <- bxp + stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", 
      tip.length = 0.01, 
      hide.ns = TRUE,
      y.position = max(df_filtered_cell_type$percentage)+5,        # Starting position for the first label
      step.increase = 0.1    # Smaller step increase to bring labels closer together
    )
  }

}

combined_plot <- wrap_plots(plot_list, ncol = 4)

pdf("figures/Fig3C.pdf", height = 9, width = 12)
combined_plot
dev.off()
