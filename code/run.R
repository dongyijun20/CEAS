library(Seurat)
library(dplyr)
merge <- readRDS("result/merged.rds")
### cell cycle analysis (some cells may be dividing)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge <- CellCycleScoring(merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merge <- ScaleData(merge, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(merge))
merge <- RunPCA(merge, npcs = 50)
merge <- merge %>% RunHarmony("sample_ident", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution=0.5) %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

saveRDS(merge,file = "result/merged.rds")

SampleColor <- readRDS("color/SampleColor.rds")
GroupColor <- readRDS("color/GroupColor.rds")
LocationColor <- readRDS("color/LocationColor.rds")
PatientColor <- readRDS("color/PatientColor.rds")
pdf(file="plot_new/merged.pdf")
DimPlot(merge, label = T)
DimPlot(merge, group.by = "sample_ident", cols = SampleColor, shuffle = T)
DimPlot(merge, group.by = "group", cols = GroupColor, shuffle = T)
DimPlot(merge, group.by = "patient_ident", cols = PatientColor, shuffle = T)
DimPlot(merge, group.by = "location", cols = LocationColor, shuffle = T)
ggplot(as.data.frame(table(merge$sample_ident)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(merge$group)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = GroupColor)+theme_bw()+xlab("Group")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
ggplot(as.data.frame(table(merge$location)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = LocationColor)+theme_bw()+xlab("Location")+ylab("Number of cells")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
DimPlot(merge, group.by = "celltype_bped_main")
DimPlot(merge, group.by = "celltype_hpca_main")
dev.off()