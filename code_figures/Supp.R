## Fig2A. plot IBD reference---------
scIBD_subset <- readRDS("reference/scIBD/scIBD_subset4use.rds")

pdf("figures/Sup_Fig2A.pdf", width = 5, height = 4)
DimPlot(scIBD_subset, label = T, raster=TRUE)
DimPlot(scIBD_subset, reduction = "umap", group.by = "major_cluster", label = TRUE, label.size = 3, raster=TRUE)
DimPlot(scIBD_subset, reduction = "umap", group.by = "disease", label = TRUE, label.size = 3, raster=TRUE)
dev.off()

pdf("figures/Sup_Fig2B.pdf", width = 16, height = 6)
DimPlot(scIBD_subset, reduction = "umap", group.by = "minor_cluster", raster=TRUE)
dev.off()