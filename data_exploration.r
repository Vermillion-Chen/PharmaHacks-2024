# install.packages("devtools")
library(devtools)
# install.packages("Seurat")
library(Seurat)
# devtools::install_github("SydneyBioX/scFeatures")
library(scFeatures)

data = readRDS("./stephenson.rds")

head(data)

