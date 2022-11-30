# load data
sparse_RNA = readRDS('/home/asmauger/biostat625final/sparse_RNA.rds')
sparse_prot = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')

all.equal(dim(sparse_RNA), c(22085, 70988))
all.equal(dim(sparse_prot), c(140, 70988))

# extract gene symbol from rownames
# rownames look like this : "ENSG00000081913-PHLPP1" (ensembl id - gene symbol)
library(stringr)
names = rownames(sparse_RNA)
hyphen_symbols = str_extract(names, "_.*")
gene_sym = str_sub(hyphen_symbols, 2)

# change rownames to just gene symbols
rownames(sparse_RNA) = gene_sym


# create seurat object
library(Seurat)
metadata = read.csv('/home/asmauger/biostat625final/metadata.csv', header=T)
rownames(metadata) = metadata[,1]
seuratobj = CreateSeuratObject(counts = sparse_RNA, meta.data=metadata)
prot_assay = CreateAssayObject(counts=sparse_prot)
seuratobj[['prot']] = prot_assay

# remove original data
rm(sparse_RNA) # access with seuratobj@assays$RNA@counts
rm(sparse_prot) # access with seuratobj@assays$prot@counts


## QC ??? not sure what we want to do here
VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# nFeature_RNA is number of unique genes expressed in a given cell
# nCount_RNA is number of transcripts in a given cell
rna_counts = seuratobj@assays$RNA@counts
all.equal(colSums(rna_counts), seuratobj$nCount_RNA)
rowSums(rna_counts!=0) # number of cells for which a gene has nonzero expression

gene_sym[grep('MT', gene_sym, perl=T)]


seuratobj@assays$prot@counts

## Normalization
seuratobj = NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(seuratobj) = 'prot'
# not sure if this is best way to normalize protein
# check non-human proteins
seuratobj = NormalizeData(seuratobj, normalization.method = 'CLR', margin = 2)

## Variable features
DefaultAssay(seuratobj) = 'RNA'
seuratobj = FindVariableFeatures(seuratobj)

## Scale data
seuratobj =ScaleData(seuratobj)
DefaultAssay(seuratobj) = 'prot'
seuratobj = ScaleData(seuratobj)

## PCA
seuratobj = RunPCA(seuratobj, assay='RNA')
DefaultAssay(seuratobj) = 'prot'
VariableFeatures(seuratobj) = rownames(seuratobj) # use all proteins
# should the Mouse and rat Igs be included??
seuratobj = RunPCA(seuratobj, assay='prot', reduction.name='protpca')

## UMAP
ElbowPlot(seuratobj)
ElbowPlot(seuratobj, reduction = 'protpca')
seuratobj = RunUMAP(seuratobj, dims=1:20)
seuratobj = RunUMAP(seuratobj, dims=1:10, reduction='protpca', reduction.name='protumap')

p1 = DimPlot(seuratobj, reduction='umap', group.by='day', shuffle = T)
p2 = DimPlot(seuratobj, reduction='umap', group.by='cell_type', shuffle=T)
p3 = DimPlot(seuratobj, reduction='protumap', group.by='day', shuffle=T)
p4 = DimPlot(seuratobj, reduction='protumap', group.by='cell_type', shuffle=T)
library(patchwork)
p1 + p2 + p3 + p4
