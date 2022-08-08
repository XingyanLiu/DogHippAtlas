library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

###########################[ function definition ]#########################


LoadCounts = function(dir_data, name,
                      dn_raw_mtx = sprintf("%s/%s_mtx", dir_data, name),
                      fn_raw_rds = sprintf("%s/matrix.rds", dn_raw_mtx),
                      gene.column = 1
){
  
  if(file.exists(fn_raw_rds)){
    message(sprintf("Loading data from %s", fn_raw_rds))
    cnt = readRDS(fn_raw_rds)
  }else{
    
    message(sprintf("Loading data from %s", dn_raw_mtx))
    cnt = Seurat::Read10X(dn_raw_mtx, gene.column = gene.column)
    
    message("backup .rds file...")
    saveRDS(cnt, fn_raw_rds)
  }
  return(cnt)
}

topKmarkers = function(markertb, ntop=5, uniq=T){
  markertb = markertb %>% group_by(cluster) %>%
    top_n(ntop, -p_val_adj) %>%
    top_n(ntop, avg_logFC)  # unnecessary
  if (uniq) {
    markers = as.character(unique(markertb$gene))
    return(markers)
  }
  return(markertb)
}

WrapperDotPlot = function(obj, 
                          genes_dot,
                          groupby=NULL,
                          gene_labs = NULL,
                          tt = NULL,
                          color_range = c("lightgrey", "mediumblue"),
                          x_rotate=90,
                          transpose = F,
                          dir_fig = NULL,
                          sname="temp",
                          w=NULL, h=NULL,
                          wscale=1., hscale=1.,
                          ...){
  message("Dot Feather plot...")
  # genes_dot = unique(markers_selected$gene)
  
  n_grps = nlevels(Idents(obj))
  
  plt_dot = DotPlot(obj, features = genes_dot, group.by = groupby,
                    cols = color_range, ...) + 
    theme(axis.text.x=element_text(angle=x_rotate, hjust=1)) #+ coord_flip()
  if (transpose){
    plt_dot = plt_dot + coord_flip()
    if(is.null(w)){w = min(5 + 0.16 * n_grps, 49.5)}
    if(is.null(h)){h = min(4 + 0.12 *length(genes_dot), 49.5)}
  }else{
    if(is.null(w)){w = min(6 + 0.15 *length(genes_dot), 49.5)}
    if(is.null(h)){h = min(4.8 + 0.16 * n_grps, 49.5)}
  }
  w = w * wscale
  h = h * hscale
  if (! is.null(gene_labs)){
    plt_dot = plt_dot + scale_x_discrete(labels = rev(gene_labs))
  }
  if (! is.null(tt)){
    plt_dot = plt_dot + ggtitle(tt)
  }
  if(!is.null(dir_fig)){
    filename = sprintf("%s/dot_%s.pdf", dir_fig, sname)
    ggsave(filename = filename, 
           plot = plt_dot, width = w,
           height = h)
    message(sprintf("figure saved into: %s size=(%.2f, %.2f)", filename, w, h))
  }
  return(plt_dot)
}

################################################################################
# "ENSCAFG00000030472 -> NRXN3
DATADIR = "E:/Users/xyliu/data003/dog"
datadir = file.path(DATADIR, "formal-hipp")

resdir = sprintf("%s/DE/%s", DATADIR, "20210902")
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)

# cnts = LoadCounts(datadir, "afterQC")
# metadata = read.csv(file.path(datadir, 'metadata.csv'), row.names = 1)
#
# #===== 0 make Seurat object & normalize =====
# obj = CreateSeuratObject(cnts, meta.data = metadata)
# obj$cluster = obj$leiden_anno
# obj = NormalizeData(obj, scale.factor = median(obj$nCount_RNA))
# str(obj)
# saveRDS(obj, file.path(datadir, "seurat_obj.rds"))

#===== read data =====
obj = readRDS(file.path(datadir, "seurat_obj.rds"))

# reorder cluster-levels
cluster_order = c(
  "0 Glutamatergic neurons",
  "3 Glutamatergic neurons",
  "5 Glutamatergic neurons",
  "6 Glutamatergic neurons",
  "7 Glutamatergic neurons",
  "10 Glutamatergic neurons",   
  "11 Glutamatergic neurons", 
  "13 Glutamatergic neurons",
  "14 Glutamatergic neurons",
  "15 Glutamatergic neurons",
  "17 Glutamatergic neurons",
  "20 Glutamatergic neurons",
  "4 GABAergic neurons",
  "8 GABAergic neurons",
  "12 GABAergic neurons",
  "21 GABAergic neurons",
  "25 Cajal-Retzius",
  "1 Astrocytes",
  "24 Microglia",
  "9 Oligodendrocyte precursor cells",
  "16 Unknown",
  "2 Myelinating oligodendrocytes",
  "18 Myelinating oligodendrocytes",
  "19 Non-neuron",
  "22 Non-neuron",
  "23 Endothelial cells"
)
obj$cluster = factor(obj$cluster, levels = cluster_order)
levels(obj$cluster)

# 2: DE test (ttest, )
groupby_de = 'cluster'  # 'leiden'
Idents(obj) <- groupby_de
tests = c('MAST', 'wilcox', 't')  #, 'roc')
test_use = tests[1]

for(test_use in tests){
  print(test_use)
  markers_all = FindAllMarkers(obj, only.pos = T, 
                               # logfc.threshold = 0.25,
                               test.use = test_use)
  
  ptag = sprintf("%s(%s)", groupby_de, test_use)
  fn_mktb = sprintf("%s/DEGtable_%s.csv", resdir, ptag)
  
  write.csv(markers_all, fn_mktb)
  print(fn_mktb)
  
  # top 10, 20, 30 markers
  for (ntop in c(50, 30, 20, 10)){
    topmk = topKmarkers(markers_all, ntop = ntop)
    write.table(topmk, sprintf("%s/top-%s-%s.csv", resdir, ntop, ptag), 
                row.names = F, col.names = F, quote = F)
  }
  
  ## TODO: handle cluster orders
  genes_dot = topKmarkers(markers_all, ntop = 5)
  WrapperDotPlot(obj, genes_dot, groupby=groupby_de, dir_fig = figdir, sname=ptag, wscale = 1.1)
  # DotPlot(obj, features = genes_dot, group.by = groupby_de) + 
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1), validate = TRUE)
  
}












