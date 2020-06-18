library(gplots)

readRDS("SNCombzCl3.rds")
alveoli <- c("AT2","CapEC","AT1","AM","FB1","FB2","matFB","MyoFB","Peri")
Alveolus <- subset(SNCombzCl3, idents = alveoli)
saveRDS(Alveolus, file = "Alveolus.rds")

bb<-read.csv("all.pairs.csv")
cluster.averages3 <- AverageExpression(Alveolus)
cluster.averages4 <- as.data.frame(AverageExpression(object = Alveolus))
##change cluster.averages3 to  count matrix. Better to use average expression value for each group.###
Alveolus.markers.2 <- FindAllMarkers(Alveolus, only.pos = TRUE, min.pct = 0.1,thresh.use = 0.1)
###change the biopsy.markers to the differential genes in your TRAP data####
biopsy.markers4<-Alveolus.markers.2[Alveolus.markers.2$avg_logFC>0.75,] ###this is to filter the de genes by log fold change

bb$Ligand.ApprovedSymbol <- tolower(bb$Ligand.ApprovedSymbol)
bb$Ligand.ApprovedSymbol <- firstup(bb$Ligand.ApprovedSymbol)

Alv_genes3<-intersect(biopsy.markers4$gene, bb$Receptor.ApprovedSymbol)
AlvReceptor3<-cluster.averages4[Alv_genes3,]
ckd_genes2
AlvReceptor3<-as.matrix(AlvReceptor3)

tiff('~/desktop/receptor_pseudotime333.tiff', units="in", width=4, height=12, res=300)
heatmap.2 (
  scale(AlvLigand2),
  col=colorRampPalette(c("yellow","lemonchiffon1","darkgreen"))(dim(AlvLigand2)[1]*dim(AlvLigand2)[2]),
  density ="none",
  Colv = F,
  Rowv = F,
  key = TRUE,
  trace ="none",
  scale="row",
  cexRow = 1,
  cexCol = 1,
  margins = c(4,10),
  keysize=1.5,
  density.info=c("histogram"),
  key.title = T,
  srtCol=45
)
dev.off()
?heatmap.2




AlveolusAvgLigands <- AverageExpression(Alveolus, return.seurat = T, features = ligands.Mouse)
AlveolusAvgReceptors <- AverageExpression(Alveolus, return.seurat = T, features = c.receptors.Mouse2)

DoHeatmap(AlveolusAvgLigands, features = c("Bdnf","Pdgfa","Wnt3a","Vegfa","Il18","Lpl","Fgf1"))
DoHeatmap(AlveolusAvgReceptors, features = c("Ngfr","Sirpa","Kdr","Lrp1","Fgf1r","Ntrk2","Hhip")) 
write.csv(AlvReceptor3, "AlvReceptors3.csv")
write.csv(AlvLigand2, "AlvLigands2.csv")

write.csv(Alveolus.markers.2, "AlvMarkers2.csv")

ALvLigand3 <- read.csv("AlvLigands low thresh.csv")
AlvReceptor4 <- read.csv("AlvReceptors low thresh.csv")


rownames(AlvReceptor4) <- AlvReceptor4[,1]
AlvReceptor4[,1] <- NULL



