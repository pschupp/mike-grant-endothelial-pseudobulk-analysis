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
- [x] Perform pseudobulk
- [x] Do FM of pseudobulk
- [x] Perform enrichments of networks
- [x] Perform DE in snRNA-seq

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
TtestOut = future_apply(dat, 1, function(datX) t.test(as.numeric(datX[which(clusters$clust == 'Endothelial cells')+1]), as.numeric(datX[-c(1,(which(clusters$clust == 'Endothelial cells')+1))])))
tvalOut = lapply(TtestOut, function(x) x$statistic)
# summary(unlist(tvalOut))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01357 0.15771 0.27415 0.47228 0.99998
out = data.frame(genes = genes$x, tvalue = unlist(tvalOut))
write.csv(out, row.names = FALSE, file = 'snrnaseq-endothelial-tvalues.csv')
```

# Log: making figure

steps:

- [x] identify module with highest enrichment
- [x] read in kme and ME
- [x] get color scheme
- [x] create correlation plot and ME plot, over each other
- [x] create ME va actual abundance plot
- [x] create t-value v kME plot
- [ ] align all plots into single figure

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
```

## Creating all plots
moving forward with whitesmoke module from "Bicor-None_signum0.147_minSize8_merge_ME_0.8_14761"

# Bicor-None_signum0.278_minSize8_merge_ME_0.8_10000,navajowhite3
```{.r}
# colorscheme
# main pink color is e465a1
# use set1 from colorbrewer
WD = "~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/10k_genes/"
network = "pseudobulk_snrna-seq_Modules/Bicor-None_signum0.278_minSize8_merge_ME_0.8_10000"
WD = "~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/"
network = "pseudobulk_snrna-seq_Modules/Bicor-None_signum0.147_minSize8_merge_ME_0.8_14761"
library('data.table')
library('ggplot2')
library('RColorBrewer')
library('future')
library('future.apply')
# load pseudobulk expression matrix
setwd(WD)
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
setwd(paste0(WD, network))
# enrich = fread(list.files(path = '.', pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE))
modSelec = 'whitesmoke'
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = fread(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE))
corRankGenes = kmeTab$Gene[order(kmeTab[,which(colnames(kmeTab) == paste0('kME', modSelec)), with = FALSE], decreasing = TRUE)]
corRankGenes = corRankGenes[-grep('^ENSG|^LINC|^ERCC|-AS1|-AS2', corRankGenes)]
corExpr = rnaDataframe[match(corRankGenes, rnaDataframe$Gene),]
corSelec = corExpr[seq(1,15),]
corSelec = data.frame(Gene = corSelec[,1], t(apply(log2(corSelec[,-1]), 1, scale)))
meSelec = data.frame(Sample = meTab$Sample, meTab[, which(colnames(meTab) == modSelec), with = F])
corPlot = reshape2::melt(corSelec, id.var = 'Gene')
colnames(corPlot) = c('Gene', 'Sample', 'Expr')
meanExpr = aggregate(corPlot$Expr, by = list(corPlot$Gene), mean)
meanExpr = meanExpr[order(meanExpr$x, decreasing = TRUE),]
corPlot$Gene = factor(corPlot$Gene, levels = meanExpr$Group.1)
corPlot$Sample = factor(corPlot$Sample, levels = unique(corPlot$Sample))
# cor between ME and actual abundance
setwd(WD)
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
clusters = data.frame(fread(list.files(path = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/', pattern = 'cluster', full.names = TRUE)))
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4', 'Endothelial cells', 'Clone 4', 'Oligodendrocytes', 'Oligodendrocytes', 'Microglia', 'Neurons')
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
abund$Cell.name = clusters$clust[match(clusters$name, abund$Cell.name)]
abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.name), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
abundPlot = data.frame(section = meSelec[,1], abundance = abundAgg[,7], me = meSelec[,2])
# abundPlot = reshape2::melt(abundPlot[order(apply(abundPlot[-1], 1, sum), decreasing = FALSE), ])
# abundPlot$section = factor(abundPlot$section, levels = unique(abundPlot$section))
# load differential expression values
library('data.table')
expr = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-expr.csv')
genes = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-gene.csv', header = T)
clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
dat = data.frame(Genes = genes[,-1], expr)
colnames(dat)[1] = 'Gene'
clusters = data.frame(clusters)
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
clusters$clust = gsub('Clone 4:1', 'Clone 4', clusters$clust)
clusters$clust = gsub('Clone 4:2', 'Clone 4', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:1', 'Oligodendrocytes', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:2', 'Oligodendrocytes', clusters$clust)
dat = dat[-grep('^ENSG|^ERCC|^LINC', dat$Gene),]
gSum = apply(dat[,-1], 1, median)
dat = dat[order(gSum, decreasing = TRUE),]
dat = data.frame(dat[1:10000,])
TtestOut = future_apply(dat[,-1] , 1, function(exprX)                                                               
    t.test(as.numeric(exprX[which(clusters$clust == 'Endothelial cells')+1]),                                     
    as.numeric(exprX[-c(1,(which(clusters$clust == 'Endothelial cells')+1))])))
tvalOut = lapply(TtestOut, function(x) x$statistic)  
deTval = data.frame(genes = dat$Gene, tvalue = unlist(tvalOut))
deTval =read.csv('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrnaseq-endothelial-tvalues.csv')
difExplot = data.frame(gene = kmeTab$Gene, kme = (kmeTab$kMEwhitesmoke), tval = (deTval$tvalue[match(kmeTab$Gene, deTval$genes)]))
# enrich = fread(list.files(path = '.', pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE))
setwd(WD)
WD = "/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/pseudobulk_snrna-seq_Modules"
modSelec = 'navajowhite3'
kmeTabbig = lapply(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE),fread)
out=lapply(kmeTabbig, function(kmeTab) apply(kmeTab[,seq(9,ncol(kmeTab),2), with=F], 2, function(x) cor(x, deTval$tvalue[match(kmeTab$Gene, deTval$genes)], use='pairwise.complete.obs')))
lapply(out,max)
# all plots made below
setwd(WD)
# create theme
fig1Theme = function(){
    theme_bw() +
    theme(
#       axis.line = element_line(colour = "black"),
# 		legend.title = element_text(size=30, family='NimbusSan'),
 		axis.text.x = element_text(size=15, color='black',  family='NimbusSan'), # , margin=margin(t=10)),
 		axis.text.y = element_text(size=15, color='black', family='NimbusSan'), # , margin=margin(r=10)),
        axis.title.y = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.x = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=0, b=0, l=0)),
        plot.title = element_text(size=30,face="bold", hjust=.5, family='NimbusSan'), # margin=margin(t=-20, b=10)),
# 		plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=10, b=10)),
# 		axis.line.x = element_line(size=3),
# 		axis.line.y = element_line(size=3),
# 		plot.margin = unit(c(4, 2, 1, 2), "lines"),
# 		legend.key.size=unit(1.3, 'cm'),
# 		legend.text=element_text(size=30, family='NimbusSan')
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')
}
# correlation plot
pdf('correlation_genes_plot.pdf')
ggplot(corPlot, aes(x = Sample, y = Expr, group = Gene)) +
    geom_line(aes(color = Gene)) +
    scale_color_discrete(colorRampPalette( brewer.pal(9,"Set1") )(15)) + # revisit TODO
    scale_x_discrete(breaks = c()) +
    scale_y_continuous(breaks = c()) +
    labs(x = '', y = 'Expression level', title = 'Pseudobulk gene\ncoexpression module') +
    guides(color=guide_legend(title="", nrow = 3, byrow = TRUE)) +
    fig1Theme()
dev.off()
# ME plot
meSelec$Sample = factor(meSelec$Sample, levels = meSelec$Sample)
pdf('ME_plot.pdf', height=2)
ggplot(meSelec, aes(x = Sample , y = lightcyan)) +
    geom_bar(stat = 'identity',fill ='#e465a1', color = 'black') +
    labs(x = 'Pseudobulk samples', y = 'AU', title = 'Module eigengene (PC1)') +
    scale_y_continuous(breaks = c(-.2, 0,.2), limits=c(-0.2,0.2)) +
    scale_x_discrete(breaks = c()) +
    fig1Theme()
dev.off()
# correlation of abundance v ME plot
cor(abundPlot$me, abundPlot$abundance)
pdf('abundance_correlation.pdf')
ggplot(abundPlot, aes(x = me, y = abundance)) +
    geom_point(color = 'black', fill = '#e465a1', size = 3, shape = 21) +
    labs(x = 'Predicted abundance in pseudobulk samples\n(module eigengene)', 
        y = 'Actual abundance in pseudobulk samples (%)', 
        title = 'Malignant cell abundance') +
    scale_y_continuous(breaks = c(0, 3,6), limits = c(0, 6)) +
    scale_x_continuous(breaks = c(-.25, 0,.25), limits = c(-.25, .25)) +
    fig1Theme()
dev.off()
# correlation between t-value and kME
pdf('difex_correlation.pdf')
ggplot(difExplot[sample(seq(1,nrow(difExplot)), 6000),], aes(x = tval, y = kme)) +
    geom_point(color = 'black', fill= '#e465a1', size = 3, shape = 21) +
    labs(x = 'Single-nucleus DE (t-values)',
        y = 'Pseudobulk kME',
        title = 'Malignant cell expression') +
#    scale_x_continuous(breaks = c(-20, 0, 20), limits = c(-20, 20)) +
#    scale_y_continuous(breaks = c(-.5, 0, 0.5, 1), limits = c(-.5, 1)) +
    fig1Theme()
dev.off()
```

```{.r}
## Read in GSEA functions, e.g.:
renv::install('bioc::GSEABase')
renv::install('bioc::limma')
renv::install('bioc::future.apply')
library('flashClust')
library('parallel')
plan(multisession, workers = 20)
source("/home/patrick/code/oldham-lab/GSEA/GSEAfxsV3.lint.R")
#source("/home/patrick/code/oldham-lab/GSEA/GSEAfxsV3.R")
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/10k_genes/pseudobulk_snrna-seq_Modules'
setwd(WD)
print("Set working directory")
## To run enrichment analysis for our gene sets in all networks:
MyGSHGloop(kmecut1="topmodposbc",exclude="none",pvalcut1=NULL)
setwd(WD)
MyGSHGloop(kmecut1="topmodposfdr",exclude="none",pvalcut1=NULL)
setwd(WD)
```


# Log 09/22

Main issues encountered with pseudobulk + FM approach are the inability to successfully recapture the relevant cell-type as modules. 

Potential issues that I will address in edit:

- too many genes, remove ENSG and lowly expressed genes
- not enough variability, boost variablity
- too few cells. will have to implement bootstrap sampling, i.e. with replacement in future version of pseudobulk code

## Pseudobulk

```{.r}
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
expr = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-expr.csv')
genes = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-gene.csv', header = T)
clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
dat = data.frame(Genes = genes[,-1], expr)
colnames(dat)[1] = 'Gene'
clusters = data.frame(clusters)
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
clusters$clust = gsub('Clone 4:1', 'Clone 4', clusters$clust)
clusters$clust = gsub('Clone 4:2', 'Clone 4', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:1', 'Oligodendrocytes', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:2', 'Oligodendrocytes', clusters$clust)
dat = dat[-grep('^ENSG|^ERCC|^LINC', dat$Gene),]
gSum = apply(dat[,-1], 1, median)
dat = dat[order(gSum, decreasing = TRUE),]
dat = data.frame(dat[1:10000,])
makeSyntheticDatasets(
    expr = dat,
    sampleindex = c(2:ncol(dat)),
    cell.info = clusters,
    cell.name = 1,
    cell.type = 3,
    cell.frac = NULL,
    pcnt.cells = 30,
    pcnt.var = 0,
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
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
plan(multisession, workers = 20)
# source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
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
    signumvec = c(.999, .99,0.95, 0.90, 0.80),
    minsizevec = c(5, 8, 10, 12, 15, 20),
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

## Check to run after FM

```{.r}
library('data.table')
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
clusters = fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/snrna-seq-cluster.csv', header = T)
colnames(clusters) = c('name', 'cluster')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clusters$clust = names(clustConv)[match(clusters$clust, clustConv)]
clusters$clust = gsub('Clone 4:1', 'Clone 4', clusters$clust)
clusters$clust = gsub('Clone 4:2', 'Clone 4', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:1', 'Oligodendrocytes', clusters$clust)
clusters$clust = gsub('Oligodendrocytes:2', 'Oligodendrocytes', clusters$clust)
abund$Cell.name = clusters$clust[match(clusters$name, abund$Cell.name)]

abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.name), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
lapply(seq_along(meTab), function(i) apply(cor(abundAgg, meTab[[i]][,-1]), 1, max, na.rm = TRUE))

# winner
# "./Bicor-None_signum0.147_minSize8_merge_ME_0.8_14761/Module_eigengenes_12-58-30.csv", whitesmoke
```

## identify module with highest enrichment

```{.r}
library('data.table')
WD = '~/investigations/sc_glioma/darmanis/GBM_data_and_metadata/'
expr = fread(paste0(WD, 'GBM_normalized_gene_counts.csv'))
clusters = fread(paste0(WD, 'GBM_metadata.csv'))
colnames(expr)[1] = 'Gene'
clusters = clusters[match(colnames(expr)[-1], clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
clusters$V1 = make.names(clusters$V1)
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/10k_genes')
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
lapply(seq_along(meTab), function(i) apply(cor(abundAgg, meTab[[i]][,-1]), 1, max, na.rm = TRUE))
cor(abundAgg, meTab[[12]][,-1])
# winner
# Bicor-None_signum0.278_minSize8_merge_ME_0.8_10000,navajowhite3
```

```{.r}
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/10k_genes/'
setwd(WD)
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
enrich = lapply(list.files(path = 'pseudobulk_snrna-seq_Modules', pattern = 'FDR.*csv', recursive = TRUE, full.names = TRUE), fread)
lapply(enrich, function(x) min(x[grep('kelley.*endothelial', x$SetName, ignore.case = TRUE), -seq(1,7)]))
#coral1 in 
modSelec = colnames(enrich)[which.min(enrich[which(enrich$SetID== 'MOSET7058'), seq(8, ncol(enrich)), with = FALSE])]
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
setwd(paste0(WD, enrichPath[netSelec]))
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
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
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
apply(cor(abundAgg, meTab[[1]]), 1, max, na.rm = TRUE)
```

## Conclusion

- using the pcnt.var was the issue. it resulted in everything being covaries with variation in the amount of cells being included.
- using fewer genes results in worse correlations between the ME and abundance

# Log 09/27

Using same methodology as above but now with Darmanis et al. dataset.

## Pseudobulk

```{.r}
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis')
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
WD = '~/investigations/sc_glioma/darmanis/GBM_data_and_metadata/'
expr = fread(paste0(WD, 'GBM_normalized_gene_counts.csv'))
clusters = fread(paste0(WD, 'GBM_metadata.csv'))
colnames(expr)[1] = 'Gene'
clusters = clusters[match(colnames(expr)[-1], clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
clusters$V1 = make.names(clusters$V1)
clusters = clusters[clusters$Location == 'Tumor',]
clusters = clusters[!(clusters$Selection == 'Unpanned'), ]
expr = expr[,c(1, match(clusters$V1, colnames(expr)[-1])+1)]
# zeros = apply(expr[,-1], 1, function(x) length(which(x == 0)))
# expr = expr[zeros<3482,]
makeSyntheticDatasets(
    expr = expr,
    sampleindex = c(2:ncol(expr)),
    cell.info = clusters,
    cell.name = 1,
    cell.type = 3,
    cell.frac = NULL,
    pcnt.cells = 30,
    pcnt.var = 0,
    no.samples = 100,
    no.datasets = 1
)
```

## FindModules

```{.r}
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis')
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
plan(multisession, workers = 20)
# source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
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
    signumvec = c(.999, .99,0.95, 0.90, 0.80),
    minsizevec = c(5, 8, 10, 12, 15, 20),
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

Network is generated under `/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis/pseudobulk_snrna-seq_Modules`.


## Check after FM

```{.r}
library('data.table')
WD = '~/investigations/sc_glioma/darmanis/GBM_data_and_metadata/'
expr = fread(paste0(WD, 'GBM_normalized_gene_counts.csv'))
clusters = fread(paste0(WD, 'GBM_metadata.csv'))
colnames(expr)[1] = 'Gene'
clusters = clusters[match(colnames(expr)[-1], clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
clusters$V1 = make.names(clusters$V1)
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis')
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
lapply(seq_along(meTab), function(i) apply(cor(abundAgg, meTab[[i]][,-1]), 1, max, na.rm = TRUE))

# winner
# Bicor-None_signum0.108_minSize10_merge_ME_0.8_21171, cyan
```


## Plot generation

steps:

- [x] identify module with highest enrichment
- [x] read in kme and ME
- [x] get color scheme
- [x] create correlation plot and ME plot, over each other
- [x] create ME va actual abundance plot
- [x] create t-value v kME plot
- [ ] align all plots into single figure

## Creating all plots
moving forward with `cyan` module from `Bicor-None_signum0.108_minSize10_merge_ME_0.8_21171`

```{.r}
# colorscheme
# main pink color is e465a1
# use set1 from colorbrewer
WD = "~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis/"
network = "pseudobulk_snrna-seq_Modules/Bicor-None_signum0.108_minSize10_merge_ME_0.8_21171/"
modSelec = 'cyan'
library('data.table')
library('ggplot2')
library('RColorBrewer')
# load pseudobulk expression matrix
setwd(WD)
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
setwd(paste0(WD, network))
# enrich = fread(list.files(path = '.', pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE))
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = fread(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE))
corRankGenes = kmeTab$Gene[order(kmeTab[,which(colnames(kmeTab) == paste0('kME', modSelec)), with = FALSE], decreasing = TRUE)]
# corRankGenes = corRankGenes[-grep('^ENSG|^LINC|^ERCC|-AS1|-AS2', corRankGenes)]
corExpr = rnaDataframe[match(corRankGenes, rnaDataframe$Gene),]
corSelec = corExpr[seq(1,15),]
corSelec = data.frame( Gene = corSelec[,1], log2(corSelec[,-1]))
meSelec = data.frame(Sample = meTab$Sample, meTab[, which(colnames(meTab) == modSelec), with = F])
corPlot = reshape2::melt(corSelec, id.var = 'Gene')
colnames(corPlot) = c('Gene', 'Sample', 'Expr')
meanExpr = aggregate(corPlot$Expr, by = list(corPlot$Gene), mean)
meanExpr = meanExpr[order(meanExpr$x, decreasing = TRUE),]
corPlot$Gene = factor(corPlot$Gene, levels = meanExpr$Group.1)
corPlot$Sample = factor(corPlot$Sample, levels = unique(corPlot$Sample))
# cor between ME and actual abundance
setwd(WD)
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))

setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis')
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
library('future')
library('future.apply')
WD = '~/investigations/sc_glioma/darmanis/GBM_data_and_metadata/'
expr = fread(paste0(WD, 'GBM_normalized_gene_counts.csv'))
clusters = fread(paste0(WD, 'GBM_metadata.csv'))
colnames(expr)[1] = 'Gene'
clusters = clusters[match(colnames(expr)[-1], clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
clusters$V1 = make.names(clusters$V1)
clusters = clusters[clusters$Location == 'Tumor',]
clusters = clusters[!(clusters$Selection == 'Unpanned'), ]
expr = expr[,c(1, match(clusters$V1, colnames(expr)[-1])+1)]

abundAgg = aggregate(abund[,-c(1)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
abundPlot = data.frame(section = meSelec[,1], abundance = abundAgg[,2], me = meSelec[,2])
# abundPlot = reshape2::melt(abundPlot[order(apply(abundPlot[-1], 1, sum), decreasing = FALSE), ])
# abundPlot$section = factor(abundPlot$section, levels = unique(abundPlot$section))
# load differential expression values
TtestOut = future_apply(expr[,-1] , 1, function(exprX)                                                               
    t.test(as.numeric(exprX[which(clusters$Selection == 'Endothelial(BSC)')+1]),                                     
    as.numeric(exprX[-c(1,(which(clusters$Selection == 'Endothelial(BSC)')+1))])))
tvalOut = lapply(TtestOut, function(x) x$statistic)  
deTval = data.frame(genes = expr$Gene, tvalue = unlist(tvalOut))
difExplot = data.frame(gene = kmeTab$Gene, kme = (kmeTab$kMEwhitesmoke), tval = (deTval$tvalue[match(kmeTab$Gene, deTval$genes)]))
# all plots made below
setwd(WD)
# create theme
fig1Theme = function(){
    theme_bw() +
    theme(
#       axis.line = element_line(colour = "black"),
# 		legend.title = element_text(size=30, family='NimbusSan'),
 		axis.text.x = element_text(size=15, color='black',  family='NimbusSan'), # , margin=margin(t=10)),
 		axis.text.y = element_text(size=15, color='black', family='NimbusSan'), # , margin=margin(r=10)),
        axis.title.y = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.x = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=0, b=0, l=0)),
        plot.title = element_text(size=30,face="bold", hjust=.5, family='NimbusSan'), # margin=margin(t=-20, b=10)),
# 		plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=10, b=10)),
# 		axis.line.x = element_line(size=3),
# 		axis.line.y = element_line(size=3),
# 		plot.margin = unit(c(4, 2, 1, 2), "lines"),
# 		legend.key.size=unit(1.3, 'cm'),
# 		legend.text=element_text(size=30, family='NimbusSan')
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')
}
# correlation plot
pdf('correlation_genes_plot.pdf')
ggplot(corPlot, aes(x = Sample, y = Expr, group = Gene)) +
    geom_line(aes(color = Gene)) +
    scale_color_discrete(colorRampPalette( brewer.pal(9,"Set1") )(15)) + # revisit TODO
    scale_x_discrete(breaks = c()) +
    scale_y_continuous(breaks = c()) +
    labs(x = '', y = 'Expression level', title = 'Pseudobulk gene\ncoexpression module') +
    guides(color=guide_legend(title="", nrow = 3, byrow = TRUE)) +
    fig1Theme()
dev.off()
# ME plot
pdf('ME_plot.pdf', height=2)
ggplot(meSelec, aes(x = Sample , y = cyan)) +
    geom_bar(stat = 'identity',fill ='#e465a1', color = 'black') +
    labs(x = 'Pseudobulk samples', y = 'AU', title = 'Module eigengene (PC1)') +
    scale_y_continuous(breaks = c(-.2, 0,.2)) +
    scale_x_discrete(breaks = c()) +
    fig1Theme()
dev.off()
# correlation of abundance v ME plot
pdf('abundance_correlation.pdf')
ggplot(abundPlot, aes(x = me, y = abundance)) +
    geom_point(color = 'black', fill = '#e465a1', size = 3, shape = 21) +
    labs(x = 'Predicted abundance in pseudobulk samples\n(module eigengene)', 
        y = 'Actual abundance in pseudobulk samples (%)', 
        title = 'Malignant cell abundance') +
#     scale_y_continuous(breaks = c(0, 4,8), limits = c(0, 8)) +
#     scale_x_continuous(breaks = c(-.3, 0,.3), limits = c(-.3, .3)) +
    fig1Theme()
dev.off()
# correlation between t-value and kME
pdf('difex_correlation.pdf')
ggplot(difExplot, aes(x = tval, y = kme)) +
    geom_point(color = 'black', fill= '#e465a1', size = 3, shape = 21) +
    labs(x = 'Single-nucleus DE (t-values)',
        y = 'Pseudobulk kME',
        title = 'Malignant cell expression') +
#     scale_x_continuous(breaks = c(-20, 0, 20), limits = c(-20, 20)) +
#     scale_y_continuous(breaks = c(-.5, 0, 0.5, 1), limits = c(-.5, 1)) +
    fig1Theme()
dev.off()
```

# Log 09/28

Using same methodology as above but now with Kreigstein dataset.

## Pseudobulk

```{.r}
WD = '~/investigations/sc_glioma/kreigstein/'
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
setwd(WD) 
expr = fread(paste0(WD, 'exprMatrix.tsv'))
expr = data.frame(Gene = expr[,1], apply(expr[,-1], 2, function(x) x/sum(x))*1E6)
clusters = fread(paste0(WD, 'meta.tsv'))
commonCells = intersect(colnames(expr)[-1], clusters$V1)
clusters = clusters[match(commonCells, clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
expr = expr[,c(1, match(commonCells, colnames(expr)[-1])+1)]
# zeros = apply(expr[,-1], 1, function(x) length(which(x == 0)))
# expr = expr[zeros<3482,]
makeSyntheticDatasets(
    expr = expr,
    sampleindex = c(2:ncol(expr)),
    cell.info = clusters,
    cell.name = 1,
    cell.type = 16,
    cell.frac = NULL,
    pcnt.cells = 30,
    pcnt.var = 0,
    no.samples = 100,
    no.datasets = 1
)
```

## FindModules

```{.r}
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/kriegstein/')
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
plan(multisession, workers = 20)
# source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
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
    signumvec = c(.999, .99,0.95, 0.90, 0.80),
    minsizevec = c(5, 8, 10, 12, 15, 20),
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

Network is generated under `/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/darmanis/pseudobulk_snrna-seq_Modules`.


## Check after FM

```{.r}
WD = '~/investigations/sc_glioma/kreigstein/'
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
setwd(WD) 
expr = fread(paste0(WD, 'exprMatrix.tsv'))
expr = data.frame(Gene = expr[,1], apply(expr[,-1], 2, function(x) x/sum(x))*1E6)
clusters = fread(paste0(WD, 'meta.tsv'))
commonCells = intersect(colnames(expr)[-1], clusters$V1)
clusters = clusters[match(commonCells, clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
expr = expr[,c(1, match(commonCells, colnames(expr)[-1])+1)]
setwd('~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/kriegstein')
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
x=lapply(seq_along(meTab), function(i) apply(cor(abundAgg, meTab[[i]][,-1]), 1, max, na.rm = TRUE))
lapply(x, function(y) y[7])
cor(meTab[[5]], abundAgg[,7])
# winner
# Bicor-None_signum0.0864_minSize8_merge_ME_0.8_17402, violet
```


## Plot generation

steps:

- [x] identify module with highest enrichment
- [x] read in kme and ME
- [x] get color scheme
- [x] create correlation plot and ME plot, over each other
- [x] create ME va actual abundance plot
- [x] create t-value v kME plot
- [ ] align all plots into single figure

## Creating all plots
moving forward with `cyan` module from `Bicor-None_signum0.108_minSize10_merge_ME_0.8_21171`

```{.r}
# colorscheme
# main pink color is e465a1
# use set1 from colorbrewer
WD = "~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/kriegstein/"
network = "pseudobulk_snrna-seq_Modules/Bicor-None_signum0.0864_minSize8_merge_ME_0.8_17402"
modSelec = 'violet'
library('data.table')
library('ggplot2')
library('RColorBrewer')
# load pseudobulk expression matrix
setwd(WD)
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
setwd(paste0(WD, network))
# enrich = fread(list.files(path = '.', pattern = 'MY_SETS.*FDR.*csv', recursive = TRUE, full.names = TRUE))
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = fread(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE))
corRankGenes = kmeTab$Gene[order(kmeTab[,which(colnames(kmeTab) == paste0('kME', modSelec)), with = FALSE], decreasing = TRUE)]
# corRankGenes = corRankGenes[-grep('^ENSG|^LINC|^ERCC|-AS1|-AS2', corRankGenes)]
corExpr = rnaDataframe[match(corRankGenes, rnaDataframe$Gene),]
corSelec = corExpr[seq(1,15),]
corSelec = data.frame( Gene = corSelec[,1], log2(corSelec[,-1]))
meSelec = data.frame(Sample = meTab$Sample, meTab[, which(colnames(meTab) == modSelec), with = F])
corPlot = reshape2::melt(corSelec, id.var = 'Gene')
colnames(corPlot) = c('Gene', 'Sample', 'Expr')
meanExpr = aggregate(corPlot$Expr, by = list(corPlot$Gene), mean)
meanExpr = meanExpr[order(meanExpr$x, decreasing = TRUE),]
corPlot$Gene = factor(corPlot$Gene, levels = meanExpr$Group.1)
corPlot$Sample = factor(corPlot$Sample, levels = unique(corPlot$Sample))
# cor between ME and actual abundance
setwd(WD)
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))

WD = '~/investigations/sc_glioma/kreigstein/'
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
setwd(WD) 
expr = fread(paste0(WD, 'exprMatrix.tsv'))
expr = data.frame(Gene = expr[,1], apply(expr[,-1], 2, function(x) x/sum(x))*1E6)
clusters = fread(paste0(WD, 'meta.tsv'))
commonCells = intersect(colnames(expr)[-1], clusters$V1)
clusters = clusters[match(commonCells, clusters$V1),]
expr = data.frame(expr)
clusters = data.frame(clusters)
expr = expr[,c(1, match(commonCells, colnames(expr)[-1])+1)]

abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
abundPlot = data.frame(section = meSelec[,1], abundance = abundAgg[,7], me = meSelec[,2])
# abundPlot = reshape2::melt(abundPlot[order(apply(abundPlot[-1], 1, sum), decreasing = FALSE), ])
# abundPlot$section = factor(abundPlot$section, levels = unique(abundPlot$section))
# load differential expression values
library('future')
library('future.apply')
options(future.globals.maxSize= +Inf)
TtestOut = future_apply(expr[,-1] , 1, function(exprX)                                                               
    t.test(as.numeric(exprX[which(clusters$Cell.Type.Assignment == 'Endothelial')+1]),                                     
    as.numeric(exprX[-c(1,(which(clusters$Cell.Type.Assignment == 'Endothelial')+1))])))
tvalOut = lapply(TtestOut, function(x) x$statistic)  
deTval = data.frame(genes = expr$gene, tvalue = unlist(tvalOut))
difExplot = data.frame(gene = kmeTab$Gene, kme = (kmeTab$kMEviolet), tval = (deTval$tvalue[match(kmeTab$Gene, deTval$genes)]))
# all plots made below
WD = "~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/kriegstein/"
setwd(WD)
# create theme
fig1Theme = function(){
    theme_bw() +
    theme(
#       axis.line = element_line(colour = "black"),
# 		legend.title = element_text(size=30, family='NimbusSan'),
 		axis.text.x = element_text(size=15, color='black',  family='NimbusSan'), # , margin=margin(t=10)),
 		axis.text.y = element_text(size=15, color='black', family='NimbusSan'), # , margin=margin(r=10)),
        axis.title.y = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.x = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=0, b=0, l=0)),
        plot.title = element_text(size=30,face="bold", hjust=.5, family='NimbusSan'), # margin=margin(t=-20, b=10)),
# 		plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=10, b=10)),
# 		axis.line.x = element_line(size=3),
# 		axis.line.y = element_line(size=3),
# 		plot.margin = unit(c(4, 2, 1, 2), "lines"),
# 		legend.key.size=unit(1.3, 'cm'),
# 		legend.text=element_text(size=30, family='NimbusSan')
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')
}
# correlation plot
pdf('correlation_genes_plot.pdf')
ggplot(corPlot, aes(x = Sample, y = Expr, group = Gene)) +
    geom_line(aes(color = Gene)) +
    scale_color_discrete(colorRampPalette( brewer.pal(9,"Set1") )(15)) + # revisit TODO
    scale_x_discrete(breaks = c()) +
    scale_y_continuous(breaks = c()) +
    labs(x = '', y = 'Expression level', title = 'Pseudobulk gene\ncoexpression module') +
    guides(color=guide_legend(title="", nrow = 3, byrow = TRUE)) +
    fig1Theme()
dev.off()
# ME plot
pdf('ME_plot.pdf', height=2)
ggplot(meSelec, aes(x = Sample , y = violet)) +
    geom_bar(stat = 'identity',fill ='#e465a1', color = 'black') +
    labs(x = 'Pseudobulk samples', y = 'AU', title = 'Module eigengene (PC1)') +
    scale_y_continuous(breaks = c(-.2, 0,.2)) +
    scale_x_discrete(breaks = c()) +
    fig1Theme()
dev.off()
# correlation of abundance v ME plot
pdf('abundance_correlation.pdf')
ggplot(abundPlot, aes(x = me, y = abundance)) +
    geom_point(color = 'black', fill = '#e465a1', size = 3, shape = 21) +
    labs(x = 'Predicted abundance in pseudobulk samples\n(module eigengene)', 
        y = 'Actual abundance in pseudobulk samples (%)', 
        title = 'Malignant cell abundance') +
#     scale_y_continuous(breaks = c(0, 4,8), limits = c(0, 8)) +
#     scale_x_continuous(breaks = c(-.3, 0,.3), limits = c(-.3, .3)) +
    fig1Theme()
dev.off()
# correlation between t-value and kME
pdf('difex_correlation.pdf')
ggplot(difExplot, aes(x = tval, y = kme)) +
    geom_point(color = 'black', fill= '#e465a1', size = 3, shape = 21) +
    labs(x = 'Single-nucleus DE (t-values)',
        y = 'Pseudobulk kME',
        title = 'Malignant cell expression') +
#     scale_x_continuous(breaks = c(-20, 0, 20), limits = c(-20, 20)) +
#     scale_y_continuous(breaks = c(-.5, 0, 0.5, 1), limits = c(-.5, 1)) +
    fig1Theme()
dev.off()
```

# Log 10/02

Final attempt to create best possible figure 1 is to use Li et al. from DLPFC. It has more than 500 endothelial cells so is a good candidate for this analysis. Will try with counts and normalized counts.

## Pseudobulk

```{.r}
WD = '/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/'
library('data.table')
expr = fread(paste0(WD, 'matrix.mtx'), skip = 2 )
colnames(expr) = c('gene', 'cell', 'count')
exprW = dcast(expr, gene ~ cell, value.var = 'count')
genes = fread(paste0(WD, 'genes.tsv'), header = FALSE)
genes_hugo = fread(paste0(WD, 'genes_mapped.tsv'), header = FALSE)
genes = genes_hugo$V2[match(genes$V2, genes_hugo$V1)]
barcs = fread(paste0(WD, 'barcodes.tsv'), header = FALSE)
barcs_anno = fread(paste0(WD, 'author_barcode_annotations.csv'), header = TRUE)
colnames(exprW) = c('Gene', make.names(barcs_anno$Cell_Class, unique = TRUE))
exprW$Gene = genes
exprW[is.na(exprW)] = 0
exprW = data.frame(exprW)
barcs_anno = data.frame(cell_name = colnames(exprW)[-1], barcs_anno)
# optional normalization
exprW = data.frame(Gene = exprW[,1], apply(exprW[,-1], 2, function(x) x/sum(x))*1E6)
# 
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/normalized'
setwd(WD) 
makeSyntheticDatasets(
    expr = exprW,
    sampleindex = c(2:ncol(exprW)),
    cell.info = barcs_anno,
    cell.name = 1,
    cell.type = 4,
    cell.frac = NULL,
    pcnt.cells = 10,
    pcnt.var = 0,
    no.samples = 100,
    no.datasets = 1
)
```

## FindModules

```{.r}
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/normalized'
setwd(WD)
library('data.table')
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
plan(multisession, workers = 20)
# source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
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
    signumvec = c(.999, .99,0.95, 0.90, 0.80),
    minsizevec = c(5, 8, 10, 12, 15, 20),
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

## Check after FM

```{.r}
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/'
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
library('data.table')
setwd(WD) 

abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
meTab = lapply(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE), fread)
x=lapply(seq_along(meTab), function(i) apply(cor(abundAgg, meTab[[i]][,-1]), 1, max, na.rm = TRUE))
lapply(x, function(y) y[7])
cor(meTab[[21]], abundAgg[,2])
# winner - no norm
# Bicor-None_signum0.293_minSize5_merge_ME_0.8_23476, brown
```


```{.r}
## Read in GSEA functions, e.g.:
library('flashClust')
library('parallel')
source("/home/patrick/code/oldham-lab/GSEA/enrichment_functions.R")
mods = unlist(unique(kmeTab[,6])[-1])
kmeTab[,6][is.na(kmeTab[,6])]="sdf"
out = enrichment_man(lapply(mods, function(x) kmeTab$Gene[kmeTab[,6]==x]), kmeTab$Gene, '/home/shared/genesets/genesets_slim')
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/'
network = 'pseudobulk_snrna-seq_Modules/'
setwd(paste0(WD, network))
MyGSHGloop(kmecut1="topmodposfdr",exclude="none",pvalcut1=NULL)
```

## Plot generation

steps:

- [x] identify module with highest enrichment
- [x] read in kme and ME
- [x] get color scheme
- [x] create correlation plot and ME plot, over each other
- [x] create ME va actual abundance plot
- [x] create t-value v kME plot
- [ ] align all plots into single figure

## Creating all plots

```{.r}
# colorscheme
# main pink color is e465a1
# use set1 from colorbrewer
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/'
network = 'pseudobulk_snrna-seq_Modules/Bicor-None_signum0.293_minSize5_merge_ME_0.8_23476/'
modSelec = 'brown'
library('data.table')
library('ggplot2')
library('RColorBrewer')
# load pseudobulk expression matrix
setwd(WD)
rnaDataframe = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'samples_[0-9]', full.names = TRUE)))
colnames(rnaDataframe)[1] = 'Gene'
rnaDataframe$Gene = make.unique(rnaDataframe$Gene)
setwd(paste0(WD, network))
kmeTab = fread(list.files(path = '.', pattern = 'kME_table.*csv', recursive = TRUE, full.names = TRUE))
meTab = fread(list.files(path = '.', pattern = 'Module_eigengenes.*csv', recursive = TRUE, full.names = TRUE))
corRankGenes = kmeTab$Gene[order(kmeTab[,which(colnames(kmeTab) == paste0('kME', modSelec)), with = FALSE], decreasing = TRUE)]
corExpr = rnaDataframe[match(corRankGenes, rnaDataframe$Gene),]
corSelec = corExpr[seq(1,15),]
corSelec = data.frame( Gene = corSelec[,1], log2(corSelec[,-1]))
meSelec = data.frame(Sample = meTab$Sample, meTab[, which(colnames(meTab) == modSelec), with = F])
corPlot = reshape2::melt(corSelec, id.var = 'Gene')
colnames(corPlot) = c('Gene', 'Sample', 'Expr')
meanExpr = aggregate(corPlot$Expr, by = list(corPlot$Gene), mean)
meanExpr = meanExpr[order(meanExpr$x, decreasing = TRUE),]
corPlot$Gene = factor(corPlot$Gene, levels = meanExpr$Group.1)
corPlot$Sample = factor(corPlot$Sample, levels = unique(corPlot$Sample))
# cor between ME and actual abundance
setwd(WD)
abund = data.frame(fread(list.files(path = 'SyntheticDatasets', pattern = 'legend', full.names = TRUE)))


WD = '/home/shared/scsn.expr_data/human_expr/postnatal/li_2018/'
library('data.table')
expr = fread(paste0(WD, 'matrix.mtx'), skip = 2 )
colnames(expr) = c('gene', 'cell', 'count')
exprW = dcast(expr, gene ~ cell, value.var = 'count')
genes = fread(paste0(WD, 'genes.tsv'), header = FALSE)
genes_hugo = fread(paste0(WD, 'genes_mapped.tsv'), header = FALSE)
genes = genes_hugo$V2[match(genes$V2, genes_hugo$V1)]
barcs = fread(paste0(WD, 'barcodes.tsv'), header = FALSE)
barcs_anno = fread(paste0(WD, 'author_barcode_annotations.csv'), header = TRUE)
colnames(exprW) = c('Gene', make.names(barcs_anno$Cell_Class, unique = TRUE))
exprW$Gene = genes
exprW[is.na(exprW)] = 0
exprW = data.frame(exprW)
barcs_anno = data.frame(cell_name = colnames(exprW)[-1], barcs_anno)
expr = exprW

abundAgg = aggregate(abund[,-c(1,2)], by = list(abund$Cell.type), sum)
tName = abundAgg[,1]
abundAgg = t(data.frame(abundAgg[,1], apply(abundAgg[,-1], 2, function(x) 100*x/sum(x)))[,-1])
colnames(abundAgg) = tName
abundPlot = data.frame(section = meSelec[,1], abundance = abundAgg[,2], me = meSelec[,2])
# abundPlot = reshape2::melt(abundPlot[order(apply(abundPlot[-1], 1, sum), decreasing = FALSE), ])
# abundPlot$section = factor(abundPlot$section, levels = unique(abundPlot$section))
# load differential expression values
library('future')
library('future.apply')
options(future.globals.maxSize= +Inf)
TtestOut = future_apply(expr[,-1] , 1, function(exprX)                                                               
    t.test(as.numeric(exprX[which(barcs_anno$Cell_Class == 'Endo')+1]),                                     
    as.numeric(exprX[-c(1,(which(barcs_anno$Cell_Class == 'Endo')+1))])))
tvalOut = lapply(TtestOut, function(x) x$statistic)  
deTval = data.frame(genes = expr$Gene, tvalue = unlist(tvalOut))
difExplot = data.frame(gene = kmeTab$Gene, kme = (kmeTab$kMEbrown), tval = (deTval$tvalue[match(kmeTab$Gene, deTval$genes)]))
# all plots made below
WD = '~/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/li/'
setwd(WD)
# create theme
fig1Theme = function(){
    theme_bw() +
    theme(
#       axis.line = element_line(colour = "black"),
# 		legend.title = element_text(size=30, family='NimbusSan'),
 		axis.text.x = element_text(size=15, color='black',  family='NimbusSan'), # , margin=margin(t=10)),
 		axis.text.y = element_text(size=15, color='black', family='NimbusSan'), # , margin=margin(r=10)),
        axis.title.y = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.x = element_text(size=20, face='bold', family='NimbusSan'), #, margin=margin(t=0, r=0, b=0, l=0)),
        plot.title = element_text(size=30,face="bold", hjust=.5, family='NimbusSan'), # margin=margin(t=-20, b=10)),
# 		plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=10, b=10)),
# 		axis.line.x = element_line(size=3),
# 		axis.line.y = element_line(size=3),
# 		plot.margin = unit(c(4, 2, 1, 2), "lines"),
# 		legend.key.size=unit(1.3, 'cm'),
# 		legend.text=element_text(size=30, family='NimbusSan')
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')
}
# correlation plot
pdf('correlation_genes_plot.pdf')
ggplot(corPlot, aes(x = Sample, y = Expr, group = Gene)) +
    geom_line(aes(color = Gene)) +
    scale_color_discrete(colorRampPalette( brewer.pal(9,"Set1") )(15)) + # revisit TODO
    scale_x_discrete(breaks = c()) +
    scale_y_continuous(breaks = c()) +
    labs(x = '', y = 'Expression level', title = 'Pseudobulk gene\ncoexpression module') +
    guides(color=guide_legend(title="", nrow = 3, byrow = TRUE)) +
    fig1Theme()
dev.off()
# ME plot
pdf('ME_plot.pdf', height=2)
ggplot(meSelec, aes(x = Sample , y = brown)) +
    geom_bar(stat = 'identity',fill ='#e465a1', color = 'black') +
    labs(x = 'Pseudobulk samples', y = 'AU', title = 'Module eigengene (PC1)') +
    scale_y_continuous(breaks = c(-.2, 0,.2)) +
    scale_x_discrete(breaks = c()) +
    fig1Theme()
dev.off()
# correlation of abundance v ME plot
pdf('abundance_correlation.pdf')
ggplot(abundPlot, aes(x = me, y = abundance)) +
    geom_point(color = 'black', fill = '#e465a1', size = 3, shape = 21) +
    labs(x = 'Predicted abundance in pseudobulk samples\n(module eigengene)', 
        y = 'Actual abundance in pseudobulk samples (%)', 
        title = 'Malignant cell abundance') +
#     scale_y_continuous(breaks = c(0, 4,8), limits = c(0, 8)) +
#     scale_x_continuous(breaks = c(-.3, 0,.3), limits = c(-.3, .3)) +
    fig1Theme()
dev.off()
# correlation between t-value and kME
pdf('difex_correlation.pdf')
ggplot(difExplot, aes(x = tval, y = kme)) +
    geom_point(color = 'black', fill= '#e465a1', size = 3, shape = 21) +
    labs(x = 'Single-nucleus DE (t-values)',
        y = 'Pseudobulk kME',
        title = 'Malignant cell expression') +
#     scale_x_continuous(breaks = c(-20, 0, 20), limits = c(-20, 20)) +
#     scale_y_continuous(breaks = c(-.5, 0, 0.5, 1), limits = c(-.5, 1)) +
    fig1Theme()
dev.off()
```

