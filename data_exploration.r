install.packages("devtools")
install.packages("Seurat")
library(devtools)
devtools::install_github("SydneyBioX/scFeatures")

library(Seurat)
library(scFeatures)
getwd()
setwd("./Coding/PharmaHacks-2024") 

cts = readRDS("./seurat.rds")
cts = UpdateSeuratObject(cts)

str(cts)

# Remove all assays except for RNA
cts[["HTOsPool1"]] <- NULL
cts[["HTOsPool2"]] <- NULL
cts[["HTOsPool3"]] <- NULL
cts[["HTOsPool4"]] <- NULL
cts[["HTOsPool5"]] <- NULL
cts[["HTOsPool6"]] <- NULL

View(cts@meta.data)
head(cts, 20)
ctsrna = cts[["RNA"]]

VlnPlot(cts, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
plot1 <- FeatureScatter(cts, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(cts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

cts <- FindVariableFeatures(cts, selection.method="vst", nfeatures =2000)
top10 <- head(VariableFeatures(cts), 10)
plot3 <- VariableFeaturePlot(cts)
LabelPoints(plot=plot3, points=top10, repel=TRUE)

all.genes <- rownames(cts)
cts <- ScaleData(cts, features=all.genes)

cts <- RunPCA(cts, features=VariableFeatures(object=cts))
print(cts[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cts, dims = 1:2, reduction = "pca")
ElbowPlot(cts)

cts <- FindNeighbors(cts, dims = 1:18)
cts <- FindClusters(cts, resolution = 0.8)
DimPlot(cts, group.by = "RNA_snn_res.0.8", label = TRUE)

head(Idents(cts), 5)

cts <- RunUMAP(cts, dims = 1:18)
DimPlot(cts, reduction = "umap")
saveRDS(cts, file = "./reduced.rds")

data = GetAssayData(object = cts, assay = "RNA", layer = "data")
head(data)


red = readRDS("./reduced.rds")
red = UpdateSeuratObject(red)

str(red)

reddata <- red@assays$RNA@data

alldata <- scFeatures:::formatData(data = reddata, celltype = red$cells, sample = red$donor)

ffeatures <- run_gene_prop_celltype(
  data = alldata,
  type = "scrna",
  num_top_gene = 100,
  ncores = 1
)
ffeatures2 <- run_gene_mean(
  alldata,
  type = "scrna", num_top_gene = 150, ncores = 1
)

head(ffeatures)
outcome = subset(red, red$disease_stage)


write.csv(ffeatures, "./final.csv", row.names=TRUE)