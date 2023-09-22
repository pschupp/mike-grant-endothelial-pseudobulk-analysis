# Log to generate figure emphasizing pseudobulk relationship

## Email prompt
> For my upcoming R01 (due Oct. 5), I would like to recreate the attached figure with a focus on endothelial cells.  This requires a few things:
> 
> 1. Identifying a SC/SN glioma dataset that includes endothelial cells.  I think the easiest place to start might be Patrick’s snRNA-seq dataset from our paper?  Also open to other suggestions.
> 2. Re-running the pseudobulk analysis and confirming that the bottom row of figure panels looks as expected when focusing on endothelial cell coexpression / differential expression.
> 3. Making cosmetic changes to the figure to introduce vasculature to panel a and update all remaining panels accordingly (with Daniel’s help).

## Goal of analysis

Can you please run it on your snRNA-seq dataset (maybe you have already done this?) and share the results?  For reference, I have copied below the code I used to build the pseudobulk dataset from the Suva scRNA-seq dataset.  (LMK if you need any of these files or other info.)
Overall goal of this analysis is to generate the following series of plots for endothelial cells:

1. Pseudobulk gene coexpression module
2. Plot of actual endothelial cell abundance v pseudobulk endothelial module ME
3. Plot of single cell DE v pseudobulk endothelial module ME

Therefore we need to perform the following analyses:
[x] Perform pseudobulk
[x] Do FM of pseudobulk
[x] Perform enrichments of networks
[x] Perform DE in snRNA-seq

# Log : data generation

Begin documenting work.

## Pseudobulk

```{.r}
renv::init()
renv::install('data.table')
# restart r
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
expr = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-expr.csv')
genes = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-gene.csv', header = T)
clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clusters$clust = gsub('Clone 4:1', 'Clone 4', clusters$clust)
clusters$clust = gsub('Clone 4:2', 'Clone 4', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:1', 'Oligodendrocytes', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:2', 'Oligodendrocytes', clusters$clust)
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
dat = data.frame(Genes = genes[,-1], expr)
colnames(dat)[1] = 'Genes'
clusters = data.frame(clusters)
makeSyntheticDatasets(
    expr = dat,
    sampleindex = c(2:ncol(dat)),
    cell.info = clusters,
    cell.name = 1,
    cell.type = 3,
    cell.frac = NULL,
    pcnt.cells = 50,
    pcnt.var = 50,
    no.samples = 100,
    no.datasets = 1
)

# quick analysis
x = read.csv(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE))
y = aggregate(x[,-c(1,2)], by = list(x$Cell.type), sum)
y=data.frame(y[,1], apply(y[,-1], 2, function(x) 100*x/sum(x)))
y
```

Two files generated:
- `SyntheticDataset1_25pcntCells_50pcntVar_100samples_08-25-51.csv` - gene x sample matrix
- `SyntheticDataset1_25pcntCells_50pcntVar_100samples_legend_08-25-51.csv` - cell/nucleus x sample binary matrix

## FindModules

```{.r}
renv::install('bioc::WGCNA')
renv::install('bioc::Clust')
renv::install('bioc::svMisc')
renv::install('bioc::qvalue')
renv::install('bioc::purrr')
renv::install('bioc::HiClimR')
renv::install('bioc::future')
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
plan(multisession, workers = 20)
setwd('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/')
FindModules(
    projectname = 'pseudobulk_snrna-seq',
    expr = rnaDataframe,
    geneinfo = c(1),
    sampleindex = seq(2,ncol(rnaDataframe)),
    samplegroups = as.factor(colnames(rnaDataframe)[-1]),
    subset = NULL,
    simMat = NULL,
    saveSimMat = FALSE,
    simType = 'Bicor',
    overlapType = 'None',
    TOtype = 'signed',
    TOdenom= 'min', 
    beta = 1,
    MIestimator = 'mi.mm',
    MIdisc = 'equalfreq',
    signumType = 'rel',
    iterate = TRUE,
    signumvec = c(.999,.99,.98, 0.95, 0.90, 0.80),
    minsizevec = c(3, 5, 8, 10, 12, 15, 20),
    signum = NULL,
    minSize = NULL ,
    merge.by = 'ME',
    merge.param = 0.8,
    export.merge.comp = T,
    ZNCcut = 2,
    calcSW = FALSE,
    loadTree = FALSE,
    writeKME = TRUE,
    calcBigModStat = FALSE,
    writeModSnap = TRUE
)
```

Network is generated under `~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/pseudobulk_snrna-seq_Modules`.

## Enrichment analysis

```{.r}
## Enrichment analysis runcode.
## Read in GSEA functions, e.g.:
renv::install('bioc::GSEABase')
renv::install('bioc::limma')
renv::install('bioc::future.apply')
library('flashClust')
library('parallel')
source("/home/patrick/code/oldham-lab/GSEA/GSEAfxsV3.lint.R")
source("/home/patrick/code/oldham-lab/GSEA/GSEAfxsV3.R")
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/pseudobulk_snrna-seq_Modules'
setwd(WD)
print("Set working directory")
## To run enrichment analysis for our gene sets in all networks:
MyGSHGloop(kmecut1="topmodposbc",exclude="none",pvalcut1=NULL)
setwd(WD)
MyGSHGloop(kmecut1="topmodposfdr",exclude="none",pvalcut1=NULL)
setwd(WD)

## Read in Broad gene sets (MolSigDBv3): 
print("Reading in 'broadSets'")
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
print("Success reading in 'broadSets'")
## To run enrichment analysis for broad gene sets in all networks:
### Note that kmecut1 can equal "seed", "topmodposbc" (recommended), or "topmodposfdr".
print("Beginning loop BC")
setwd(WD)
BroadGSHGloop(kmecut1="topmodposbc",pvalcut1=NULL)
setwd(WD)
print("Beginning loop FD")
BroadGSHGloop(kmecut1="topmodposfdr",pvalcut1=NULL)
setwd(WD)
```

## DE in snRNA-seq data

```{.r}
# read in expression data
library('future.apply')
library('data.table')
expr = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-expr.csv')
genes = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-gene.csv', header = T)
clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
dat = data.frame(Genes = genes[,-1],expr)
clusters = data.frame(clusters)
testOut = future_apply(dat, 1, function(datX) wilcox.test(as.numeric(datX[which(clusters$clust == 'Endothelial cells')+1]), as.numeric(datX[-c(1,(which(clusters$clust == 'Endothelial cells')+1))])))
TtestOut = future_apply(dat, 1, function(datX) t.test(as.numeric(datX[which(clusters$clust == 'Endothelial cells')+1]), as.numeric(datX[-c(1,(which(clusters$clust == 'Endothelial cells')+1))])))
tvalOut = lapply(TtestOut, function(x) x$p.value)
wilcox.test(dat[ind], x[!ind], alternative='greater')$p.value
# summary(unlist(tvalOut))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01357 0.15771 0.27415 0.47228 0.99998
out = data.frame(genes = genes$x, tvalue = unlist(tvalOut))
write.csv(out, row.names = FALSE, file = 'snrnaseq-endothelial-tvalues.csv')
```

# Log: making figure

steps:

[x] identify module with highest enrichment
[x] read in kme and ME
[x] get color scheme
[ ] create correlation plot and ME plot, over each other
[ ] create ME va actual abundance plot
[ ] create t-value v kME plot
[ ] align all plots into single figure

## identify module with highest enrichment

```{.r}
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/pseudobulk_snrna-seq_Modules/'
library('data.table')
rnaDataframe = data.frame(fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/SyntheticDatasets/SyntheticDataset1_25pcntCells_50pcntVar_100samples_08-25-51.csv'))
rnaDataframe$x = make.unique(rnaDataframe$x)
enrich = lapply(list.files(path = WD, pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE), fread)
# MOSET7058 is top 150 genes from kelley endothelial
maxEnrich = lapply(enrich, function(x) min(x[which(x$setid == 'moset7058'), seq(8, ncol(x)), with = false]))
netSelec = which.min(unlist(maxEnrich))
enrichPath = list.files(path = WD, pattern = 'Bicor-None')
netSelec = which.min(unlist(maxEnrich))
setwd(paste0(WD, enrichPath[netSelec]))
enrich = fread(list.files(path = '.', pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE))
modSelec = colnames(enrich)[which.min(enrich[which(enrich$SetID== 'MOSET7058'), seq(8, ncol(enrich)), with = FALSE])]
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = fread(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE))
# colorscheme
# main pink color is e465a1
# use set1 from colorbrewer
corRankGenes = kmeTab$x[order(kmeTab[,which(colnames(kmeTab) == paste0('kME', modSelec)), with = FALSE], decreasing = TRUE)]
corRankGenes = corRankGenes[-grep('^ENSG|^LINC|^ERCC', corRankGenes)]
corExpr = rnaDataframe[match(corRankGenes, rnaDataframe$x),]
corSelec = corExpr[seq(1,5),]
cor(t(corSelec[,-1]))
corSelec = data.frame( Gene = corSelec[,1], log2(corSelec[,-1]))
cor(t(corSelec[,-1]))
meSelec = meTab[, which(colnames(meTab) == modSelec), with = F]
library('ggplot2')
setwd('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis')
# correlation plot
corPlot = reshape2::melt(corSelec, id.var = 'Gene')
colnames(corPlot) = c('Gene', 'Sample', 'Expr')
corPlot$Sample = factor(corPlot$Sample, levels = unique(corPlot$Sample))
pdf('correlation_genes_plot.pdf')
ggplot(corPlot, aes(x = Sample, y = Expr, group = Gene)) +
    geom_line(aes(color = Gene))
dev.off()


# ME plot
pdf('ME_plot.pdf')


# cor between ME and actual abundance
abund = read.csv("/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/SyntheticDatasets/SyntheticDataset1_25pcntCells_50pcntVar_100samples_legend_08-25-51.csv")

clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
abund$Cell.name = clusters$clust[match(clusters$name, abund$Cell.name)]

abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.name), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
```
