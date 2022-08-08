
# BiocManager::install("BioQC")
library(ggplot2)
library(dplyr)
library(stringr)

library(BioBase)
library(BioQC)

############################## functions ###############################
write_namelist = function(lst, fp, quote=F){
  write.table(lst, fp,
              row.names=F, col.names=F, quote=quote)
  message(fp)
}
load_namelist = function(fp, header = F, ...){
  message(fp)
  read.delim(fp, header = header, as.is=T, ...)[,1]
}

load_multiple_lists = function(dirname, suffix='.txt', ...) {
  multilists = list()
  fnames = list.files(dirname, pattern=suffix)
  # print(fnames)
  for (.fn in fnames) {
    .nm = gsub(suffix, '', .fn)
    multilists[[.nm]] = load_namelist(sprintf("%s/%s", dirname, .fn))
  }
  multilists
}


############################## directory ###############################
MAINDIR = "E:/Users/xyliu/data003/dog"
datadir = "E:/Users/xyliu/data003/dog/formal-hipp"
GENEDIR = file.path(MAINDIR, "formal-genes")

# datadir_subset = "E:/Users/xyliu/data003/dog/hipp_subset"
# fn = file.path(datadir, "agg_expr_cutprop0.25.csv")
# exprMat = read.csv(fn, row.names = 1)  # 4523 g &  4281 mini-cl
# exprMat[1: 5, 1:6]
# saveRDS(exprMat, file.path(datadir_subset, "agg_expr_cutprop0.25.rds"))

# mat_agg = readRDS(file.path(datadir_subset, "agg_expr_cutprop0.25.rds"))
# mat_agg[1:5, 1:7]
# write_namelist(rownames(mat_agg), sprintf("%s/genes_clprop0.25.csv", datadir_subset))

###########################################################################
# load PSG lists
all_psg_lists = load_multiple_lists(sprintf("%s/PSGs", GENEDIR))
all_psg_lists %>% lapply(length)

# expression averages & proportions
expr_avgs = read.csv(sprintf("%s/expr_avgs.csv", datadir), row.names = 1)
expr_props = read.csv(sprintf("%s/expr_props.csv", datadir), row.names = 1)

expr_props_max = rowMax(as.matrix(expr_props))  # for gene filtering

cut_prop = 0.05
genes_high_prop = rownames(expr_props[expr_props_max >= cut_prop,])
print(length(genes_high_prop))

# filter out un-detected genes; and transformed into row-indices
idx_map = seq(genes_high_prop)
names(idx_map) <- genes_high_prop

all_psg_indices = list()
for (nm in names(all_psg_lists)) {
  .genes = intersect(all_psg_lists[[nm]], genes_high_prop)
  all_psg_indices[[nm]] = as.numeric(idx_map[.genes])
}
all_psg_indices %>% lapply(length)

mat_use = as(expr_avgs, 'matrix')[genes_high_prop, ]
res = wmwTest(
  mat_use, all_psg_indices, valType = "p.greater", simplify = F)
res
matrix(p.adjust(res, "BH"), ncol=dim(mat_use)[2])


library(pheatmap)
pmat = - log10(res) 
pmat = pmat / apply(pmat, 1, max)
pheatmap(pmat, cluster_rows = F, cluster_cols = F, )


############# Entropy ##############

avg_score_vs_bg = function(scores, ids_or_names) {
  ids_or_names = intersect(ids_or_names, names(scores))
  n_fg = length(ids_or_names)
  n_bg = length(scores)
  message(sprintf("n_fg = %s, n_all = %s", n_fg, n_bg))
  sum_fg = sum(scores[ids_or_names])
  fg = sum_fg / n_fg
  bg = (sum(scores, na.rm = T) - sum_fg) / (n_bg - n_fg)
  
  c("fg_score" = fg, "bg_score" = bg)
}
calc_gene_entropy_specifity = function(mat, norm=T) {
  df_entropy = data.frame(
    gene_entropy = apply(mat, 1, entropy),
    gene_specificity = entropySpecificity(mat, norm=norm)
  )
  df_entropy
}

calc_multi_score_vs_bg = function(scores, gene_lists, tags = c('PSG', 'background')) {
  res = list()
  for (tag in names(gene_lists)) {
    .genes = gene_lists[[tag]]
    res[[tag]] = avg_score_vs_bg(scores, .genes)
  }
  bind_cols(type=tags, res)
}


##=================
library(Seurat)

obj = readRDS(file.path(datadir, "seurat_obj.rds"))
str(obj@meta.data)

resdir = file.path(datadir, "entropy")
## mat_use1 = mat_use1 %>% as("Matrix") %>% as("dgCMatrix")  # test

#### major / minor class Avgs
avg_by = c("major_class", "leiden_anno", NULL)[1]
if (is.null(avg_by)) {
  #### single-cell expression matrix
  mat_base = obj[['RNA']]@data
  rownames(mat_base) = rownames(obj)
} else {
  Idents(obj) <- avg_by
  mat_base = AverageExpression(obj, assays = "RNA")[[1]]
  
}

####### compute gene scores (entropy and specificity) #############
df_entropy = calc_gene_entropy_specifity(mat_base)

summary(df_entropy)
write.csv(df_entropy, sprintf("%s/df_entropy-%s.csv", resdir, avg_by))

#===== select a proper background ====
# cut_prop1 = 0.2
# score_type = c("gene_entropy", "gene_specificity")[2]
resdf = data.frame()
for (cut_prop1 in c(0, 0.1, 0.2, 0.25)) {
  bg_genes = rownames(expr_props[expr_props_max >= cut_prop1,])
  
  # 1
  score_type = "gene_entropy"
  gene_scores = df_entropy[bg_genes, score_type]
  names(gene_scores) = bg_genes
  
  .resdf = calc_multi_score_vs_bg(
    gene_scores, all_psg_lists, 
    tags = paste(score_type, c('PSG', 'bg'), cut_prop1, sep='-')
  )
  resdf = bind_rows(resdf, .resdf)
  
  # 2
  score_type = "gene_specificity"
  gene_scores = df_entropy[bg_genes, score_type]
  names(gene_scores) = bg_genes
  .resdf = calc_multi_score_vs_bg(
    gene_scores, all_psg_lists, 
    tags = paste(score_type, c('PSG', 'bg'), cut_prop1, sep='-')
  )
  resdf = bind_rows(resdf, .resdf)
  
}

resdf

write.csv(resdf, sprintf("%s/psg_entropy_specificity-%s.csv", resdir, avg_by), row.names = F)
# hist(gene_entropy)


################################[ dot-plot ]###############################
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
# obj$cluster
obj$cluster = factor(obj$cluster, levels = cluster_order)
levels(obj$cluster)

# dot-plot
figdir = file.path(resdir, 'figs')
print(names(all_psg_lists))

nm_psg = names(all_psg_lists)[2]
genes_dot = intersect(all_psg_lists[[nm_psg]], rownames(obj))
length(genes_dot)
dot_groupby = "cluster"
WrapperDotPlot(obj, genes_dot, groupby=dot_groupby, 
               dir_fig = figdir, 
               sname=paste(dot_groupby, nm_psg, sep='-'), 
               hscale = 1.1,
               wscale = 1.25)


# Seurat::VlnPlot(obj, features = genes_dot, group.by = dot_groupby)
###########################################################


library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function for violin-plots
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          x_rotate = 45,
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=x_rotate, hjust=1, vjust = 1), 
          axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


vln_groupby = 'cluster'
Idents(obj) <- vln_groupby

nm_psg = names(all_psg_lists)[2]
genes_vln = intersect(all_psg_lists[[nm_psg]], rownames(obj))

p_vln = StackedVlnPlot(obj, genes_vln)
# p_vln
filename = sprintf("%s/vln_%s.pdf", figdir, paste(vln_groupby, nm_psg, sep='-'))
w = min(5 + 0.16 * nlevels(Idents(obj)), 49.5)
h = min(4 + 0.16 *length(genes_vln), 49.5)
ggsave(filename = filename, plot = p_vln, 
       width = w,
       height = h
       )


#######################################
fp_deg = sprintf("%s/DEGtable_cluster(MAST).csv", "E:/Users/xyliu/data003/dog/DE/20210902")
marker_df = read.csv(fp_deg, as.is=c("X", "gene"))
str(marker_df)
genes_vln = topKmarkers(marker_df, ntop = 3)

p_vln = StackedVlnPlot(obj, genes_vln)
# p_vln
filename = sprintf("%s/vln_%s.pdf", figdir, paste(vln_groupby, "top3", sep='-'))
w = min(5 + 0.16 * nlevels(Idents(obj)), 49.5)
h = min(4 + 0.16 *length(genes_vln), 49.5)
ggsave(filename = filename, plot = p_vln, 
       width = w,
       height = h
)


