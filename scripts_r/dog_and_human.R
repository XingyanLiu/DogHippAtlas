##########################################################
library(Seurat)
# library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(harmony)

WORKDIR = 'E:/lxy_pro/003/DogHippocampus/scripts_r'
# WORKDIR = "D:/Users/xyliu/005/scripts_r"  # 192-server
setwd(WORKDIR)
source('Rutils_biomics.R')
source('Rutils_align.R')

'%!in%' <- function(x, y){!('%in%'(x,y))}

load_h5_meta = function(dn, name, as_srt=T){
  fp_h5 = sprintf("%s/%s.h5", dn, name)
  fp_meta = sprintf("%s/metadata-%s.csv", dn, name)
  
  cnts = Read10X_h5(fp_h5)
  meta = read.csv(fp_meta, row.names = 1)
  if (as_srt) {
    obj = CreateSeuratObject(cnts, meta.data = meta)
    return(obj)
  } else {
    return(list(counts=cnts, metadata=meta))
  }
  
}
load_human_hipp = function(as_srt=F){
  return(load_h5_meta(
    dn="F:/external_data/human_hipp/mats",
    name = "human_hipp-sample_each_1000-leiden", as_srt=as_srt))
}
load_dog_hipp = function(as_srt=F){
  return(load_h5_meta(
    dn="F:/003_dog/mats", 
    name = "dog_hipp-sample_each_1000-leiden", as_srt=as_srt))
}
  
datadir = "F:/003_dog"
resdir = file.path(datadir, "crosssp")
figdir = file.path(resdir, 'figs')
dir.create(figdir, recursive = T)

# ==== Load homologies ====
fp_gmap = "../resources/homo-dog2human.csv"
df_gmap = read.csv(fp_gmap)
str(df_gmap[,c(2, 4)])
df_gmap = take_only_1v1(df_gmap[,c(2, 4)])
str(df_gmap)

# ==== Load counts and metadata =====
obj1 = load_dog_hipp(as_srt = F)
obj2 = load_human_hipp(as_srt = F)
obj1
obj2
obj1$metadata$species = 'dog'
obj2$metadata$species = "human"

cnts.list = align_datasets(obj1$counts, obj2$counts, df_gmap, uni_names = T)
obj1 = CreateSeuratObject(cnts.list[[1]], meta.data = obj1$metadata)
obj2 = CreateSeuratObject(cnts.list[[2]], meta.data = obj2$metadata)

merged = merge(obj1, obj2)

merged <- integrateByHarmony(
  merged, resdir=resdir, figdir=figdir, 
  key_data = 'species', key_type='major_type',
  scale_split='species'
)

Seurat::DimPlot(merged, split.by = 'species', group.by = 'major_type')
write.csv(merged@meta.data, sprintf("%s/metadata-merged.csv", resdir))
str(merged@meta.data)

