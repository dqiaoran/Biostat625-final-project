library(rhdf5)

## Load rnaseq data ##
h5ls("/home/asmauger/biostat625final/train_cite_inputs_raw.h5")
pathinput = "/home/asmauger/biostat625final/train_cite_inputs_raw.h5"
rna_data = h5read(pathinput, "/train_cite_inputs_raw/block0_values", index=list(1:22085, 1:70988))
rna_rows = h5read(pathinput, '/train_cite_inputs_raw/axis0')
rna_cols = h5read(pathinput, 'train_cite_inputs_raw/axis1')

## Load protein data ##
h5ls("/home/asmauger/biostat625final/train_cite_targets_raw.h5")
pathinput='/home/asmauger/biostat625final/train_cite_targets_raw.h5'
prot_data = h5read(pathinput, '/train_cite_targets_raw/block0_values')
prot_rows = h5read(pathinput, '/train_cite_targets_raw/axis0')
prot_cols = h5read(pathinput, '/train_cite_targets_raw/axis1')

all.equal(rna_cols, prot_cols)


## Sparse matrix

library(Matrix)
smalldata = as(rna_data, Class='dgCMatrix')
dimnames(smalldata) = list(rna_rows, rna_cols)
smalldata_prot = as(prot_data, Class='dgCMatrix')
dimnames(smalldata_prot) = list(prot_rows, prot_cols)

saveRDS(smalldata, file='sparse_RNA.rds')
saveRDS(smalldata_prot, file='sparse_prot.rds')


## Use this to load sparse matrix:

sparse_RNA = readRDS('/home/asmauger/biostat625final/sparse_RNA.rds')
sparse_prot = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')

all.equal(dim(sparse_RNA), c(22085, 70988))
all.equal(dim(sparse_prot), c(140, 70988))

## Seurat object

library(Seurat)
metadata = read.csv('/home/asmauger/biostat625final/metadata.csv', header=T)
seuratobj = CreateSeuratObject(counts = sparse_RNA, meta.data=metadata)
prot_assay = CreateAssayObject(counts=sparse_prot)
seuratobj[['prot']] = prot_assay

Assays(seuratobj) # should have RNA and prot

DefaultAssay(mydata2) = 'prot' # set default assay like so


