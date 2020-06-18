#Exon only vs intron+exon gene detection...input is zUMIs output count matrix...


reads<-readRDS("/home/humphreys/JEFF/lung_nuc2/outs/zUMIs_output/expression/SNucLung2Z_JEFF_NUC2_NUC2_NUC2.dgecounts.rds") #Load zumis count table
#select reads (intron or intron + exon)
reads<-reads$umicount$intron$all    #introns
reads <- reads$umicount$inex$all    #introns and exons
dim(reads)
intron.cells <- colnames(reads)

length(intersect(intron.cells, all.cells))

#Gene transfer
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
aa<-data.frame(listAttributes(mouse))
results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),mart = mouse)
all1<-data.frame(as.matrix(intron.cells))
all1$gene<-rownames(all1)
annotate2<-results[which(results$ensembl_gene_id %in% all1$gene),]
all1<-all1[annotate2$ensembl_gene_id,]
all1<-cbind(all1, annotate2)
all1<-all1[!duplicated(all1$external_gene_name),]
rownames(all1)<-all1$external_gene_name
Nuc1UMIin <- all1[,-c((ncol(all1)-2):ncol(all1))]
Nuc1UMIin <-Matrix(as.matrix(Nuc2UMIin), sparse = T)
saveRDS(Nuc1UMIin, "Nuc1UMIin.rds")

#Create seurat object and check nFeature, nCount
library(Seurat)
CellUMIintrexon <- CreateSeuratObject(CellUMIinex, project = "snUMI", assay = "RNA", min.cells = 5)
mean(SCellLung1z@meta.data$nFeature_RNA )
mean(SNucLung2z@meta.data$nCount_RNA)
saveRDS(Nuc1UMIintrexon, file = "Nuc1UMIintrexon.rds")
CellUMIintrexon
CellReadsIntrexon





CellUMIin <- reads$umicount$intron$all
CellReadsIn <- reads$readcount$intron$all
CellUMIex <- reads$umicount$exon$all
CellReadsex <- reads$readcount$exon$all
CellUMIinex <- reads$umicount$inex$all
CellReadsinex <- reads$readcount$inex$all
Nuc1UMIex <- reads$umicount$exon$all
Nuc1UMIin <- reads$umicount$intron$all
Nuc1UMIinex <- reads$umicount$inex$all

dim(CellReadsex)
dim(Nuc1UMIex)
dim(Nuc1ReadsEx)
dim(Nuc1UMIinex)
dim(Nuc2UMIin)
dim(CellUMIinex)

library(Seurat)
CellUMIintrexon <- CreateSeuratObject(CellUMIinex, project = "snUMI", assay = "RNA", min.cells = 5)
mean(SCellLung1z@meta.data$nFeature_RNA )
mean(SNucLung2z@meta.data$nCount_RNA)
saveRDS(Nuc1UMIintrexon, file = "Nuc1UMIintrexon.rds")
CellUMIintrexon
CellReadsIntrexon

reads<-reads[, all.cells]
dim(reads)
