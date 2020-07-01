####Subset data for Immune, Endothelial, and Epithelial Cells


#Subset Epithelial cells

Epiz2 <- subset(SNCombzCl3, idents = c("Bas","Club","Mes","AT2","AT1","Cil"))
Epiz2 <- FindVariableFeatures(Epiz2, selection.method = "vst", nfeatures = 2000)
Epiz2 <- ScaleData(Epiz2, display.progress = F, features = rownames(Epiz2))
Epiz2 <- RunPCA(Epiz2, features = VariableFeatures(object = Epiz2))
ElbowPlot(Epiz2, ndims = 50)
Epiz2 <- FindNeighbors(Epiz2, reduction = "harmony", dims = 1:30)
Epiz2 <- FindClusters(Epiz2, reduction = "harmony", resolution = 0.7)
Epiz2 <- RunUMAP(Epiz2, reduction = "harmony", dims = 1:30,min.dist = 0.08,seed.use = 20)
DimPlot(Epiz2, label = F, label.size = 3) + NoLegend()
Epiz2.markers <- FindAllMarkers(Epiz2, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
FeaturePlot(Epiz2, c("Ngfr"), order = T)
saveRDS(Epiz2, file = "Epiz2.rds")
write.csv(Epiz2.markers, file = "Epiz2Markers.csv")
new.cluster.ids <- c("AT2","AT1","AT2","Cili","Club","Meso","Gob","AT1","Basal")
names(new.cluster.ids) <- levels(Epiz2)
Epiz2 <- RenameIdents(Epiz2, new.cluster.ids)
DotPlot(Epiz2, features = EpizMarkers)
EpizMarkers <- c("Acoxl","Adam19","Il33","Hopx","Lama3","Cdh13","Dnah9","Rp1","Atp2b2",
                 "Nwd2","Rasgrf1","Sulf1","Muc16","Wdr72","Muc5b","Kcnma1","Pcsk1","Trp63")

#Generic code to plot marker lists on dotplot
DotPlot(Endoz, features = rev(EndozMarkers), cols = c("grey", "black"),
        dot.scale = 8) + RotatedAxis() + theme(axis.text=element_text(size=8),
                                               axis.title=element_text(size=8,face="bold"))


#Subset immune cells
ImmCellz <- subset(SNCombzCl3, idents = c("AM","TC","BC","Th17 TC","TC3","Mono","IM","DC","GB"))
ImmCellz <- FindVariableFeatures(ImmCellz, selection.method = "vst", nfeatures = 2000)
ImmCellz <- ScaleData(ImmCellz, display.progress = F, features = rownames(ImmCellz))
ImmCellz <- RunPCA(ImmCellz, features = VariableFeatures(object = ImmCellz))
ElbowPlot(ImmCellz, ndims = 50)
ImmCellz <- FindNeighbors(ImmCellz, reduction = "harmony", dims = 1:40)
ImmCellz <- FindClusters(ImmCellz, reduction = "harmony", resolution = 1.0)
ImmCellz <- RunUMAP(ImmCellz, reduction = "harmony", dims = 1:40)
DimPlot(ImmCellz2, label = T, label.size = 3) + NoLegend()
ImmCellz2.markers <- FindAllMarkers(ImmCellz2, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
FeaturePlot(ImmCellz2, c("C1qc"), order = T)

write.csv(ImmCellz.markers, file = "ImmCellz2markers.csv")

saveRDS(ImmCellz, file = "ImmCellz2.rds")
ImmCellz.markers <- FindAllMarkers(ImmCellz, min.pct = 0.10, logfc.threshold = 0.1)

#Subset Endothelial cells
Endoz <- subset(SNCombzCl3, idents = c("EC1","aEC","vEC","CapEC","LEC"))
Endoz <- FindVariableFeatures(Endoz, selection.method = "vst", nfeatures = 2000)
Endoz <- ScaleData(Endoz, display.progress = F, features = rownames(Endoz))
Endoz <- RunPCA(Endoz, features = VariableFeatures(object = Endoz))
ElbowPlot(Endoz, ndims = 50)
Endoz <- FindNeighbors(Endoz, reduction = "harmony", dims = 1:25)
Endoz <- FindClusters(Endoz, reduction = "harmony", resolution = 0.8)
Endoz <- RunUMAP(Endoz, reduction = "harmony", dims = 1:25)
DimPlot(Endoz, label = T, label.size = 3) + NoLegend()
FeaturePlot(Endoz, c("Il5"), order = T)
saveRDS(Endoz, file = "Endoz.rds")

Endoz.markers <- FindAllMarkers(Endoz, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
write.csv(Endoz.markers, file = "EndozMarkers.csv")

