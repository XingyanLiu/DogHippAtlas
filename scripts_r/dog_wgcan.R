# BiocManager::install("latticeExtra")
# BiocManager::install("WGCNA")
# install.packages("msigdbr")
# install.packages("vroom")
# BiocManager::install("DOSE")
# BiocManager::install("clusterProfiler")
# BiocManager::install("magrittr")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Cf.eg.db") # Genome wide annotation for Canine (dog)

library("WGCNA")
library("dplyr")
library("RColorBrewer")

options(stringsAsFactors = FALSE)
enableWGCNAThreads(8)

MAIN_DIR = "D:/Users/xyliu/003"
setwd(MAIN_DIR)
source("RFunxWGCAN.R")

DATADIR = "E:/Users/xyliu/data003/dog"

datadir = file.path(DATADIR, "hipp_subset")
list.files(datadir)

# ====================================================================
#                     1. data input and preprocessing
# ====================================================================

### settings
tag.type = c("sc", "agg")[2]
tag.sub = c("", "_sub")[1]
tag.cut = 0.25
nm_mat = sprintf("%s_expr%s_cutprop%.2g", tag.type, tag.sub, tag.cut)
nm_mat

resdir = sprintf("%s/WGCAN/%s/%s", DATADIR, nm_mat, "0802")
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)
message(sprintf("directory for saving results: \n\t%s", resdir))


### load expression matrix
ExprMat = read.csv(sprintf("%s/%s.csv", datadir, nm_mat), row.names = 1) 
ExprMat = t(ExprMat) # transposed to sample-gene matrix
# saveRDS(ExprMat, file = sprintf("%s/ExprMat.rds", resdir))

metadata = read.csv(sprintf("%s/%s_metadata.csv", datadir, nm_mat))
rownames(metadata) = rownames(ExprMat) # [1:5]
# metadata$leiden_anno[metadata$leiden == 3] = "3 Glutamatergic neurons" 
str(metadata)

ExprMat[1:5, 1:5]
n_genes = ncol(ExprMat)
n_samples = nrow(ExprMat)

##################


#=====[ Domestication genes - TEST ]=====
DomeGenes = c("AGAP1", "ANKS1B", "CCBE1", "CUX2", 
              "ENSCAFG00000022711", "GAD2", "MBP", 
              "OPCML", "PCSK5", "ROR1")
DomeGenes %in% colnames(ExprMat)

# ====================================================================
#                     2. Network construction
# ====================================================================


# ===========2.1 parameter selection===========
# choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from=12, to=20, by=2)) # powers

# call the network topology analysis function
sft = pickSoftThreshold(ExprMat, powerVector = powers, blockSize=NULL,
                        verbose = 5)

# ====plot the results====
source("RFunxWGCAN.R")
plotSft4PowerSelection(sft = sft, fn_pdf = sprintf("%s/power_selection.pdf", figdir))


# ========[ 2a.2 network construction and module detection ]==========
# TOM: Topological Overlap Matrices
param_power = 4
param_corr = "pearson"
minModuleSize = 10
mergeCutHeight = 0.55  # 0.4 # 0.25
net = blockwiseModules(ExprMat, power=param_power,
                       corType = param_corr,
                       maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = minModuleSize,
                       reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
names(net)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)#, colorSeq = )
table(moduleColors)
# net$MEs contains the module eigengenes of the modules.

####===============####
fn_net = sprintf("%s/net_pwr%d_minsz%s_mgcut%.2g.rds", 
               resdir, param_power, minModuleSize, mergeCutHeight)
saveRDS(net, file=fn_net)
# net = readRDS(fn_net)


str(moduleColors)

geneModuleLabels = as.data.frame(moduleLabels)
geneModuleLabels$gene = rownames(geneModuleLabels)
geneModuleLabels$moduleColors = moduleColors
geneModuleLabels$moduleNames = paste0("ME", moduleColors)
### needs to add gene-module membership (following code)
# write.csv(geneModuleLabels, sprintf("%s/gene_module_labels.csv", resdir))

# str(geneModuleLabels)

#==========================================================================
# plot the results
# sizeGrWindow(12, 9)
{
pdf(sprintf("%s/gene_cluster_denrograms_%.2g.pdf", figdir, cut, mergeCutHeight), 
    width = 6, height = 5)
plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()
}

### ============[ module correlations ]=================
MEs = net$MEs
ME_colors = labels2colors(as.numeric(stringr::str_remove(colnames(MEs), "ME")))
colnames(MEs) = paste0("ME", ME_colors)
MEs = orderMEs(MEs)

{
pdf(sprintf("%s/eigengene_adjacency_%.2g.pdf", figdir, mergeCutHeight),
    width = 6, height = 8)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      heatmapColors = COLORS4VALUES,
                      marDendro = c(3, 3, 2, 4),
                      marHeatmap = c(3, 4, 2, 2),
                      plotDendrograms = T)
dev.off()
}


### ==================[ module - cluster correlations (1) ]====================
### (with binary labels)

cln = 'leiden_anno'
labels_orig = metadata[[cln]]
type_ord = order(unique(metadata$leiden))
all_types = unique(labels_orig)[type_ord]
labels_binary = oneHotLabel(labels_orig, all_types = all_types)
# str(labels_binary)
type_ord_plot = read.csv(sprintf("%s/formal-hipp/celltype_ordered.csv", DATADIR), header=F)$V1

#...

### ============[ Module-GeneExpr-BinaryLabel correlations (2, better) ]====================

### MM: module membership
### GS: gene significance

geneModuleMembership = as.data.frame(cor(ExprMat, MEs, use = "p"))
rownames(geneModuleMembership)[1:5]
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), n_samples));

geneTypeSignificance = as.data.frame(cor(ExprMat, labels_binary, use = "p"))
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTypeSignificance), n_samples));
# str(GSPvalue)


moduleTypeCor = cor(geneModuleMembership, 
                    geneTypeSignificance, 
                    use="p") 

write.csv(moduleTypeCor, sprintf("%s/mod_type_corr.csv", resdir))
###[ re-order the matrix ]###

source("RFunxWGCAN.R")
# mat = clustReorderAll(moduleTypeCor)   #, distance = 'correlation')
mat = moduleTypeCor[, type_ord_plot]
row_ord = order(apply(mat, 1, which.max))
mat = mat[row_ord, ]
# col_ord = order(apply(mat, 2, which.max))
# mat = mat[, col_ord]


### remove `grey` module (it is the background module)
mat = mat[rownames(mat) != 'MEgrey', ]
write.csv(mat, sprintf("%s/mod_type_corr (mat).csv", resdir))
  
### plot heatmap
{
tt = "Module-cluster relationships"
hmapModuleTypeCor(mat, #moduleTypeCor, 
                  tt = tt,
                  width = 10, height=5.5,
                  fn_pdf = sprintf("%s/%s_%.2g.pdf", figdir, tt, mergeCutHeight))
}
table(moduleColors)

### memberships of domestication genes
# pheatmap::pheatmap(geneModuleMembership[DomeGenes, ], 
#                    cluster_rows = T, cluster_cols = T)

### save gene-module menbership matrix (gene-by-module)
write.csv(geneModuleMembership, sprintf("%s/gene_mod_membership.csv", resdir))

save(
  geneModuleMembership,
  geneTypeSignificance,
  moduleTypeCor,
  file = sprintf("%s/G_M_T_corr.RData", resdir)
)


# ======================================================================
#       arrangement of results
# ======================================================================

### adding `memvership` column to data.frame `geneModuleLabels`
geneModuleLabels$membership <- rownames(geneModuleLabels) %>% 
  sapply(function(x){
  geneModuleMembership[x, geneModuleLabels[x, "moduleNames"]]
}) 
str(geneModuleLabels)

### module summary
moduleSummary <- geneModuleLabels %>% 
  group_by(moduleNames) %>% 
  summarise(n_genes = n(), 
            median = median(membership), 
            mean = mean(membership), 
            max = max(membership))

ntop = 20
moduleTopGenes <- geneModuleLabels %>% 
  group_by(moduleNames) %>%
  top_n(membership, n = ntop) %>%
  arrange(moduleLabels, - membership)

moduleTopGenes %>% summarise(min = min(membership), max =max(membership), n = n())

cor_cut = 0.75
moduleTopGenes1 <- geneModuleLabels %>% 
  group_by(moduleNames) %>%
  subset(membership > cor_cut) %>%
  arrange(moduleLabels, - membership)
moduleTopGenes1 %>% summarise(min = min(membership), max =max(membership), n = n())

write.csv(geneModuleLabels, sprintf("%s/gene_module_labels.csv", resdir))
write.csv(geneModuleLabels[DomeGenes, ], sprintf("%s/dome_gene_module_labels.csv", resdir))

write.csv(moduleSummary, sprintf("%s/module_summary.csv", resdir), row.names = F)
write.csv(moduleTopGenes, sprintf("%s/module_genes_top%d.csv", resdir, ntop), row.names = F)
write.csv(moduleTopGenes1, sprintf("%s/module_genes_corr%.2g.csv", resdir, cor_cut), row.names = F)

lapply(unique(geneModuleLabels$moduleNames), function(x){
  subdf = subset(geneModuleLabels, moduleNames == x) %>% arrange(-membership)
  write.csv(subdf, sprintf("%s/%s_genes.csv", resdir, x), row.names = F)
})





################################################################################
###     overlaps with Top50 DEGs (for each cluster)
################################################################################

venndir = file.path(figdir, 'venn')
dir.create(venndir)

### loading DEGs
genedir = 'E:/Users/xyliu/data003/dog/formal-hipp'
deg_tb = read.csv(sprintf("%s/analyzed_marker_names.csv", genedir))
str(deg_tb)

###===================[ Define functions ]=====================
### function for computing overlap-scores
overlapScore = function(set1, set2){
  n1 = length(set1)
  n2 = length(set2)
  nn = length(intersect(set1, set2))
  s = nn / min(n1, n2)
}


### function for Venn plot
library(VennDiagram)
venn_colors = c("#f6ab6c", "#6a2c70")

# plotVennPair = function(set1, set2, )
venn.plot = venn.diagram(
  x = list(set1, set2),
  height = 2000, width = 2000,
  filename = fname,
  fill = venn_colors,
  col = venn_colors,
  margin = 0.1,
  category.names = c(name1, name2)
)



###===================================================
### preparing name-list for looking-up
TypeNames = as.list(all_types)  # Note that NO cluster-id pasted!
tmpstr = stringr::str_split_fixed(all_types, ' ', n=2)
names(TypeNames) <- paste0("X", tmpstr[, 1])

CIDs = colnames(deg_tb)
cid = CIDs[1]

for (cid in CIDs){
cname = TypeNames[[cid]]
cname

### for each cluster, extract the most correlated module
### `moduleTypeCor`, or `mat` without `grey` module
# modname = which.max(mat[, cname]) %>% names()
modname = which.max(ovs_tb[, cname]) %>% names()
modname

###===== take out module genes =======
# head(geneModuleLabels)
modgenes = subset(geneModuleLabels, moduleNames == modname)$gene
length(modgenes)

###===== take out cluster DE genes =======
ntop_deg = 50
degs_use = deg_tb[[cid]] %>% head(ntop_deg)

###===== overlap scores (ovs) =======
ovs = overlapScore(set1 = modgenes, set2 = degs_use)

print(ovs)

###======= Venn plot ========
venn_colors = c('1' = "#f6ab6c", '2' = "#6a2c70")
venn.plot = venn.diagram(
  x = list('1' = modgenes, '2' = degs_use),
  height = 2000, width = 2000,
  filename = sprintf("%s/venn_%s_%s.tif", venndir, modname, cname),
  fill = venn_colors,
  col = venn_colors,
  margin = 0.2,
  main = sprintf("overlap score: %.2f", ovs),
  category.names = c(modname, cname)
)

}

###############[ computing overlap-scores for all pairs ]###############
use_mods = rownames(mat)

sapply(colnames(mat), function(cn){
  cid = paste0("X", strsplit(cn, ' ')[[1]][1])
  
  cg = deg_tb[[cid]] %>% head(ntop_deg) # cluster DEGs
  sapply(use_mods, function(md){
    mdg = subset(geneModuleLabels, moduleNames == md)$gene # module genes
    ovs = overlapScore(mdg, cg)
    ovs
  })
}) -> ovs_tb

write.csv(ovs_tb, sprintf("%s/table_overlap_scores.csv", resdir))

hmapModuleTypeCor(ovs_tb, #moduleTypeCor,
                  tt = tt_ovs,
                  width = 10, height=5.5,
                  fn_pdf = sprintf("%s/%s_%.2g.pdf", figdir, tt_ovs, mergeCutHeight))




