library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)

merge <- readRDS("result/merged.rds")
pbmc <- merge
pbmc$batch <- pbmc$sample
pbmc$batch <- pbmc$group

table(pbmc$batch)
set.seed(10)
pbmc <- pbmc %>% 
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

pbmc <- pbmc %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 
pbmc <- pbmc %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

saveRDS(pbmc,"result/merged_rmbatch.rds")

pdf(file="plot/merged_rmbatch.pdf")
DimPlot(pbmc, label = T)
DimPlot(pbmc, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(pbmc, group.by = "group", cols = GroupColor, shuffle = T)
ggplot(as.data.frame(table(pbmc$sample)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")
DimPlot(pbmc, group.by = "celltype_bped_main", label = T)
DimPlot(pbmc, group.by = "celltype_hpca_main", label = T)
dev.off()

