library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
#library(muscat)
library(reshape2)
library(celldex)
library(clustree)
library(SingleR)
library(tidyverse)
library(RColorBrewer)

setwd("E:/Rstudio_default_working/gse")   


a1 <- Read10X(data.dir ="E:/Rstudio_default_working/gse/PBL.Pt8")

pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
#创建对象，counts 就是读取的那个矩阵，project 样本的名字，min.cells =3, min.features=200一个基因至少在三个细胞中表达才会被保留，或者一个细胞必须有200个基因才能保留
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#计算线粒体基因的百分比
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot用小提琴图来进行可视化，ncol = 3用几列数据绘制，一些死细胞和一些质量差的细胞线粒体基因比较高，某些细胞除外。nrow()行名。
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#FeatureScatter函数来可视化特征-特征之间的关系即将这些基因可视化出来，nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例。FeatureScatter表示nCount_RNA与percent.mt的相关性如为－是负相关，+为正相关
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
#表示过滤一些低质量的细胞，可以自己根据需要修改
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#数据进行标准化（标准化之前counts数是整数，标准化后变成小数），通常是取log的形式，为了消除文库深度对测序的影响，是通过总的表达值对每一个细胞的基因表达值进行进行的标准化最后在乘以比例因子（10000），之后对结果进行对数的转换。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#鉴定高度变化的基因，这些基因用于后续分析



# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)  #查看前五个高变基因
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE)
