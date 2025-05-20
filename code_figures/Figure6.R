## Fig6B: filter genes----------
endo_rm_lym <- readRDS("result/endo_rm_lym.rds")
markeres_CEAS <- FindMarkers(endo_rm_lym, group.by = "group", ident.1 = "CEAS", ident.2 = "HC", only.pos = T)
markers_cluster3 <- FindAllMarkers(endo_rm_lym, group.by = "EC_type", only.pos = T)
intersected_markers <- intersect(rownames(markeres_CEAS)[markeres_CEAS$p_val_adj < 0.01 & markeres_CEAS$avg_log2FC > 1],
                                 rownames(markers_cluster3)[markers_cluster3$cluster == "CEAS-specific" & markers_cluster3$p_val_adj < 0.01 & markers_cluster3$avg_log2FC > 1])
intersected_markers

sup_table <- markers_cluster3[markers_cluster3$cluster=="CEAS-specific"&rownames(markers_cluster3)%in%intersected_markers,] %>% 
  arrange(avg_log2FC)
write.csv(sup_table, file = "tables/Sup_Table1.csv", quote = F)

interest_genes <- c("CTNNA3","ADAMTS6","CCDC3")
interest_genes %in% intersected_markers

pdf("figures/Fig6B.pdf", height = 3, width = 9)
FeaturePlot(endo_rm_lym, interest_genes, max.cutoff = "q95", min.cutoff = "q5", cols = c("grey","coral"), ncol = 3, pt.size = 0.3)
dev.off()
