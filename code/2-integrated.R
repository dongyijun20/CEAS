library(Seurat)

CE001yz<-readRDS("result/CE001yz.rds")
CE002yz<-readRDS("result/CE002yz.rds")
CE003yz<-readRDS("result/CE003yz.rds")
CE004yz<-readRDS("result/CE004yz.rds")
CE005yz<-readRDS("result/CE005yz.rds")
CE006yz<-readRDS("result/CE006yz.rds")
CE007yz<-readRDS("result/CE007yz.rds")
CE008zc<-readRDS("result/CE008zc.rds")
CE009yz<-readRDS("result/CE009yz.rds")
CE0010hc<-readRDS("result/CE0010hc.rds")
CE0011hc<-readRDS("result/CE0011hc.rds")
CE0012hc<-readRDS("result/CE0012hc.rds")

### give label
CE001yz$location <- "Colon"
CE003yz$location <- "Colon"
CE006yz$location <- "Colon"
CE007yz$location <- "Colon"
CE008zc$location <- "Colon"
CE0010hc$location <- "Colon"
CE0012hc$location <- "Colon"

CE002yz$location <- "Ileum"
CE009yz$location <- "Ileum"
CE0011hc$location <- "Ileum"

CE004yz$location <- "Duodenum"
CE005yz$location <- "Stomach"

CE001yz$group <- "CEAS"
CE002yz$group <- "CEAS"
CE003yz$group <- "CEAS"
CE004yz$group <- "CEAS"
CE005yz$group <- "CEAS"
CE006yz$group <- "CEAS"
CE007yz$group <- "CEAS"
CE009yz$group <- "CEAS"

CE008zc$group <- "SC"

CE0010hc$group <- "HC"
CE0011hc$group <- "HC"
CE0012hc$group <- "HC"

CE001yz$sample_ident <- "Sample1"
CE002yz$sample_ident <- "Sample2"
CE003yz$sample_ident <- "Sample3"
CE004yz$sample_ident <- "Sample4"
CE005yz$sample_ident <- "Sample5"
CE006yz$sample_ident <- "Sample6"
CE007yz$sample_ident <- "Sample7"
CE008zc$sample_ident <- "Sample8"
CE009yz$sample_ident <- "Sample9"
CE0010hc$sample_ident <- "Sample10"
CE0011hc$sample_ident <- "Sample11"
CE0012hc$sample_ident <- "Sample12"

CE001yz$patient_ident <- "Patient1"
CE002yz$patient_ident <- "Patient2"
CE003yz$patient_ident <- "Patient3"
CE004yz$patient_ident <- "Patient4"
CE005yz$patient_ident <- "Patient4"
CE006yz$patient_ident <- "Patient4"
CE007yz$patient_ident <- "Patient5"
CE008zc$patient_ident <- "Patient5"
CE009yz$patient_ident <- "Patient6"
CE0010hc$patient_ident <- "Patient7"
CE0011hc$patient_ident <- "Patient8"
CE0012hc$patient_ident <- "Patient8"

# set color
library(RColorBrewer)
SampleColor <- ArchRPalettes$stallion[1:12]
names(SampleColor)<-paste0("Sample",1:12)
saveRDS(SampleColor, file = "color/SampleColor.rds")

GroupColor <- brewer.pal(3, 'Set1')
names(GroupColor) <- c("CEAS","HC","SC")
saveRDS(GroupColor, file = "color/GroupColor.rds")

LocationColor <- brewer.pal(4, 'Set2')
names(LocationColor) <- c("Colon","Ileum","Duodenum","Stomach")
saveRDS(LocationColor, file = "color/LocationColor.rds")

PatientColor <- ArchRPalettes$kelly[1:8]
names(PatientColor) <- paste0("Patient",1:8)
saveRDS(PatientColor, file = "color/PatientColor.rds")

merge <- merge(x=CE001yz,y=c(CE002yz,CE003yz,CE004yz,CE005yz,CE006yz,CE007yz,CE008zc,CE009yz,CE0010hc,CE0011hc,CE0012hc))
merge
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge)
merge <- ScaleData(merge)

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

merge$sample_ident<-as.factor(merge$sample_ident)
levels(merge$sample_ident) <- paste0("Sample",1:12)
merge$patient_ident<-as.factor(merge$patient_ident)
levels(merge$patient_ident) <- paste0("Patient",1:8)
merge$group<-as.factor(merge$group)
levels(merge$group) <- c("CEAS","HC","SC")

table(merge$seurat_clusters)
merge <- merge[,merge$seurat_clusters%in%names(table(merge$seurat_clusters))[table(merge$seurat_clusters)>500]]
merge$seurat_clusters<-as.factor(merge$seurat_clusters)
levels(merge$seurat_clusters) <- names(table(merge$seurat_clusters))[table(merge$seurat_clusters)>500]

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

table(merge$seurat_clusters, merge$celltype_bped_main)



