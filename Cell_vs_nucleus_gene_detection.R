library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

#Analysis/annotation of merged cell and nucleus object

DimPlot(NCCombz2, label = F, label.size = 3, group.by = "stim2") + NoLegend()
FeaturePlot(NCCombz2, c("Nox4"), order = T)
NCCombz2
NCCombz2.markers <- FindAllMarkers(NCCombz2, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
write.csv(NCCombz2.markers, file = "NCCombz2Markers.csv")

new.cluster.ids <- c("EC", "AT2", "AM", "FB1", "AT1", "capEC", "Club", "TC1", "Cil", "Peri",
                     "cMono", "vEC", "BC", "ncMono", "aEC", "Mes", "Th17 TC", "DC1", "Cd103+ DC",
                     "matFB", "AT2b", "IM", "GB", "TC3?", "MyoFB", "capEC2", "Nphil", "LEC",
                     "Nphil/FB", "Div", "SMC", "Cd209+/Siglech+ DC", "Div2", "Peri/EC", "AM/AT", "AM/RBC")
names(new.cluster.ids) <- levels(NCCombz2)
NCCombz2 <- RenameIdents(NCCombz2, new.cluster.ids)

#Isolation of mesenchymal cells and reclustering
NCMesenchz <- subset(NCCombz2, idents = c("FB1","Peri","matFB","MyoFB","SMC"))
NCMesenchz <- FindVariableFeatures(NCMesenchz, selection.method = "vst", nfeatures = 2000)
NCMesenchz <- ScaleData(NCMesenchz, display.progress = F, features = rownames(NCMesenchz))
NCMesenchz <- RunPCA(NCMesenchz, features = VariableFeatures(object = NCMesenchz))
ElbowPlot(NCMesenchz, ndims = 50)
NCMesenchz <- FindNeighbors(NCMesenchz, reduction = "harmony", dims = 1:40)
NCMesenchz <- FindClusters(NCMesenchz, reduction = "harmony", resolution = 0.7)
NCMesenchz <- RunUMAP(NCMesenchz, reduction = "harmony", dims = 1:40)
NCMesenchz <- RunUMAP(NCMesenchz, dims = 1:20, min.dist = 0.008,seed.use = 808)

DimPlot(NCMesenchz, label = T, label.size = 3, split.by = "stim2") + NoLegend()
FeaturePlot(NCMesenchz, c("Dcn"), order = T, split.by = "stim2")

saveRDS(NCCombz2, file = "NCCombz2.rds")

DoHeatmap(NCMesenchz, features = c("Gm38407"), group.by = "stim2")
DotPlot(NCMesenchz, features = ZincFingerTF, group.by = "stim2")
VlnPlot(NCMesenchz, features = c("Pmepa1"), pt.size = 0)


############Volcano plot etc for cell v nucleus genes ###########
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

all.data<-GetAssayData(NCMesenchz)
table(NCMesenchz@meta.data$stim2)
cell.raw.data<-as.matrix(all.data[,1:1329])
nuc.raw.data<-as.matrix(all.data[,1330:3728])

nuc.data<-data.frame(gene=rownames(nuc.raw.data),Nuc_expr=rowMeans(nuc.raw.data))
cell.data<-data.frame(gene=rownames(cell.raw.data),Cell_expr=rowMeans(cell.raw.data))

avg.expr<-merge(nuc.data,cell.data, by='gene')
rownames(avg.expr)<-avg.expr$gene

#Calc proportion of cells v nuclei expressing a gene
all.data<-data.frame(as.matrix(all.data))

######Idents -> Cell or nucleus##########
# Store the cluster identities in a new column in 'object@meta.data'
NCMesenchz$saved.idents <- Idents(object = NCMesenchz)
# Set the experimental condition as cell identity
object <- Seurat::SetAllIdent(object = object, id = "stim")
Idents(object = NCMesenchz) <- "stim2"

#####Resume volcano plot procedure...###

nuc_id<-WhichCells(object = NCMesenchz, ident = 'Nucleus')
cell_id<-WhichCells(object = NCMesenchz, ident = 'Cell')
nuc_prop <- round(
  x = apply(
    X = all.data[, nuc_id, drop = F],
    MARGIN = 1,
    FUN = function(x) {
      return(sum(x > 0) / length(x = x))
    }
  ),
  digits = 3
)
cell_prop <- round(
  x = apply(
    X = all.data[, cell_id, drop = F],
    MARGIN = 1,
    FUN = function(x) {
      return(sum(x > 0) / length(x = x))
    }
  ),
  digits = 3
)
gene_prop <- cbind(data.frame(nuc_prop), data.frame(cell_prop))
colnames(x = gene_prop) <- c("Nuc_prop","Cell_prop")


#Wilcoxon Rank Sums...Substitute Findmarkers for Seurat3....
aa <- FindMarkers(NCMesenchz, ident.1 = "Cell", logfc.threshold = 0,
            test.use = "wilcox", min.pct = 0, only.pos = FALSE,min.cells.feature = 0,
            min.cells.group = 0)

bb <- FindMarkers(NCMesenchz, ident.1 = "Nucleus", logfc.threshold = 0,
                  test.use = "wilcox", min.pct = 0, only.pos = FALSE,min.cells.feature = 0,
                  min.cells.group = 0)


p_val<-WilcoxDETest(NCMesenchz, cells.1 = nuc_id, cells.2 = cell_id) #this step may take a long while.
p_val<-na.omit(p_val)
nuc_expm <- apply(X = all.data[, nuc_id, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 1))
cell_expm <- apply(X = all.data[, cell_id, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 1))
total.diff <- (nuc_expm - cell_expm)
aa$logFC2 <- total.diff[rownames(aa)]
bb$logFC2 <- total.diff[rownames(bb)]
aa$p_val_adj2 = p.adjust(p = aa$p_val, method = "bonferroni", n = nrow(aa))

####Compute expected variation (95% CI)###

meta.data<-NCMesenchz@meta.data
meta.data$stim2 <- as.character(meta.data$stim2)
cell.ids <- which(meta.data$stim2 == "Cell")
replace.names <- sample(cell.ids, round(length(cell.ids) / 2), replace = FALSE)
meta.data$stim2[replace.names] <- "Cellsb"
meta.data$stim2 <- as.factor(meta.data$stim2)

all.data<-data.frame(as.matrix(all.data))
cell.nuc.gene <- as.data.frame(matrix(NA, nrow(all.data), 
                                      nlevels(meta.data$stim2)))
rownames(cell.nuc.gene) <- rownames(all.data)
colnames(cell.nuc.gene) <- levels(meta.data$stim2)
gene.expr <- as.data.frame(matrix(NA, nrow(all.data), 
                                  nlevels(meta.data$stim2)))
rownames(gene.expr) <- rownames(all.data)
colnames(gene.expr) <- levels(meta.data$stim2)

for (i in 1:nlevels(meta.data$stim2)) {
  ident.use <- levels(meta.data$stim2)[i]
  col.subset <- which(meta.data$stim2 == ident.use)
  cell.nuc1 <- apply(all.data[, col.subset], 1, function(x) sum(x > 0))
  cell.nuc.gene[, i] <- cell.nuc1
  expr1 <- apply(all.data[, col.subset], 1, function(x) mean(x[x > 0]))
  expr1[is.na(expr1)] <- 0
  gene.expr[, i] <- expr1
}

gene.prop <- sweep(cell.nuc.gene, 2, apply(cell.nuc.gene, 2, max), "/")
cut1 <- cut(avg.expr$Cell_expr, c(0, seq(0.01, 12, by=0.2)), include.lowest = TRUE)
prop.diff <- abs(gene.prop$Cell - gene.prop$Cellsb)
prop.exp <- data.frame(Cell_expr = c(0, seq(0.11, 11.9, 0.2)), 
                       exprq = tapply(prop.diff, cut1, quantile, 0.975))
cut2 <- cut(gene.prop$Cell, seq(0, 1, by=0.02), include.lowest = TRUE)
prop.exp2 <- data.frame(Cell_expr = seq(0.01, 0.99, 0.02), 
                        exprq = tapply(gene.prop$Cellsb, cut2, quantile, 0.975))


#######Combine data for plotting....#####
gene_prop$gene<-rownames(gene_prop)
p_val$gene<-rownames(p_val)

aa$gene<-rownames(aa)
bb$gene<-rownames(bb)

gene.expr.prop1<-merge(avg.expr, gene_prop, by='gene')
gene.expr.prop1<-merge(gene.expr.prop1, aa, by='gene')

gene.expr.prop2<-merge(avg.expr, gene_prop, by='gene')
gene.expr.prop2<-merge(gene.expr.prop2, bb, by='gene')

rownames(gene.expr.prop1)<-gene.expr.prop1$gene
rownames(gene.expr.prop2)<-gene.expr.prop2$gene



##Binned scatter plot ###
library(ggplot2)
# plot proportion of nuclei and cells for gene expression
pal.spectral <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(100)
ggplot(gene.expr.prop1, aes(x = Nuc_prop, y = Cell_prop)) +
  geom_hline(yintercept = 0) +
  stat_bin_hex(bins = 50) +
  scale_fill_gradientn(colours=pal.spectral, trans="log10", name="No. of genes") +
  xlab("Nuclei gene detection") + 
  ylab("Cells gene detection") + 
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_smooth(data = prop.exp2, aes(Cell_expr, exprq), 
              color = "grey", size = 1, se = FALSE) +
  geom_smooth(data = prop.exp2, aes(exprq, Cell_expr), 
              color = "grey", size = 1, se = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=12))

# plot differential genes and significance
ggplot(gene.expr.prop1, aes(x = logFC2, y = -log10(p_val_adj))) +
  stat_bin_hex(bins = 50) +
  geom_hline(yintercept = -log10(0.05), color = "grey", size = 0.5) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), color = "grey", size = 0.5) +
  xlim(c(-2, 2)) +
  scale_fill_gradientn(colours=pal.spectral, trans="log10", name="No. of genes") +
  xlab("Nuclei vs. Cells (logFC)") +
  ylab("Significance (-log10 P-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

 #Pull out significant genes from nucleus and cell...
table((gene.expr.prop1$logFC2)<(-0.5))
table((gene.expr.prop1$logFC2)>(0.5))

table((gene.expr.prop1$diff1)>(0.25))
table((gene.expr.prop1$diff2)>(0.25))

gene.expr.prop1$diff1 <- gene.expr.prop1$Nuc_prop - gene.expr.prop1$Cell_prop
gene.expr.prop1$diff2 <- gene.expr.prop1$Cell_prop - gene.expr.prop1$Nuc_prop
table((gene.expr.prop1$diff1)<0.2) OR ((gene.expr.prop1$diff2)<0.2)
table((gene.expr.prop1$diff2)<0.2)

nodiff <- gene.expr.prop2[gene.expr.prop2$diff1 < 0.2 & gene.expr.prop2$diff2 < 0.2, ]
nodiff <- subset(gene.expr.prop1, subset= gene.expr.prop1$diff1 < 0.2 & gene.expr.prop1$diff2 < 0.2)

write.csv(gene.expr.prop2, file = "TotalMesenchGenes.csv")
gene.expr.prop2 <- read.csv("TotalMesenchGenes.csv")

CellGenesCl <- subset(gene.expr.prop1, logFC2 < (-0.5))
NucGenesCl <- subset(gene.expr.prop1, logFC2 > 0.5)
write.csv(CellGenesCl, file = "CellGenesCl.csv")
write.csv(NucGenesCl, file = "NucGenesCl.csv")


#Select genes for heatmap

Ribosomal <- c("Rps12", "Rps14", "Rps15", "Rps24", "Rps15a", "Rps4x", "Rplp1", "Rps9", "Rplp2")

Mitochondrial <- c("mt-Co2","mt-Atp8","mt-Atp6","mt-Nd4l","mt-Co1","mt-Nd2","mt-Nd1",
                   "mt-Nd3","mt-Nd6")

HeatShock <- c("Hspa1a","Hspa1b", "Hsp90aa1","Hsp90b1","Hsp90ab1","Hspa8","Hspa5","Hspb1","Hspe1")


ZincFingerTF <- c("Zbtb16","Zfp277","Zfp280d","Zfp346","Zfp521","Zfp568","Zfpm2","Zfr","Zhx3")
ZFTFCell <- c("Zfand5","Zfp36","Zfp36l1","Zfp503")

MesGenes <- c("Myh11","Hhip", "Col14a1", "Adcy8", "Nalcn","Bmper")
 

#TF database comparisons
mouseTF <- read.csv(file = "mouseTranscriptionFactorDatabase.csv")
NCMesenchz.markers <- read.csv("NCMesenchz.markers.csv")
mouseTF <- mouseTF$Gene.name
NucGenesList <- NucGenesCl$gene
CellGenesList <- CellGenesCl$gene
MesenchGenes <- NCMesenchzCl.markers$gene
MesenchGenesCl <- subset(NCMesenchzCl.markers, subset = NCMesenchzCl.markers$avg_logFC > 0.75)
MesenchGenesCl <- MesenchGenesCl$gene
MesenchTFCl <- intersect(mouseTF, MesenchGenesCl)
NucTF <- intersect(mouseTF, NucGenesList)
CellTF <- intersect(mouseTF, CellGenesList)

write.csv(NucTF, file = "NucTFCl.csv")
write.csv(CellTF, file = "CellTFCl.csv")
write.csv(MesenchTFCl, file = "MesenchTFCl.csv")

#uMAPs etc for figures
DimPlot(NCMesenchz, label = F) + NoLegend()
DimPlot(NCMesenchz, group.by = "stim2") + NoLegend()
FeaturePlot(NCMesenchz, features = c("Pdgfra"), order = T, cols= c("gray92", "darkblue"), pt.size = 1) + NoLegend()
FeaturePlot(NCMesenchz, c("Nt5e"), order = T)

VlnPlot(NCMesenchz, c("Slit3"), pt.size = 0, group.by = "stim2") + NoLegend()
dev.off()

SNCombzCl3
DimPlot(SNCCombzCl4, label = F) + NoLegend()
SNCCombzCl4
DotPlot(NCMesenchz, features = c("Hopx","Lamp3","Foxj1","Sftpc","Sftpd","Scgb1a1",
                                   "Scgb3a2"), cols = c("grey", "black"), 
        dot.scale = 6) + RotatedAxis()


#Average clusters expression

saveRDS(NCMesenchzCl, file = "NCMesenchzClAVG.rds")

# reorder clusters
Endoz@active.ident <- factor(Endoz@active.ident, 
                            levels=c("LECs", 
                                     "aECs",
                                     "vECs", 
                                     "capEC", 
                                     "ECs"))

new.cluster.ids <- c("Bmper+ FB","Brinp1+ FB","Peri","Col14a1+ FB","MyoFB","SMC")
names(new.cluster.ids) <- levels(NCMesenchzCl)
NCMesenchzCl <- RenameIdents(NCMesenchzCl, new.cluster.ids)


MesenchDotMarkers <- c("Bmper","Fgfr4","5033421B08Rik","Upk1b","Brinp1","Nalcn","Gm38407",
                       "Pamr1","Grik4","Pdgfb","Pde5a","Notch3","Col14a1","Dcn","Lsamp","Opcml",
                       "Ebf2","Ank2","Hhip","Aspn","Mapk4","Grem2","Zfp385b","Lgr5","Rbfox1","Myh11",
                       "Myocd","Acta2","Adcy5","Stac")

MesenchDotMarkers2 <- rev(MesenchDotMarkers)


#Stacked violin plots
stress
mygene <- c("Junb","Jund","Ier3","Ier2","Btg1")

MesGenes2 <- rev(MesGenes)

plot_list = list()
for (i in 1:length(MesGenes2)) {
  p = VlnPlot(NCMesenchzCl,MesGenes2[i],pt.size = 0)
  p=p+coord_flip()+xlab("")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_blank(),plot.title = element_text(size=8, face = "bold",hjust = 0.5) , legend.position="none")+ 
    theme(plot.margin=margin(l=-0.5,unit="cm"))
  plot_list[[i]] = p
}
dev.off()
plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "hv")
