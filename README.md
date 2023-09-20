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
[ ] Perform enrichments of networks
[ ] Perform DE in snRNA-seq

## Log

Begin documenting work.

## Pseudobulk

```{.r}
renv::init()
# restart r
source('~/code/oldham-lab/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r')
renv::install('data.table')
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
makeSyntheticDatasets(
    expr = dat,
    sampleindex = c(2:ncol(dat)),
    cell.info = clusters,
    cell.name = 1,
    cell.type = 2,
    cell.frac = NULL,
    pcnt.cells = 25,
    pcnt.var = 50,
    no.samples = 100,
    no.datasets = 1
)
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
rnaDataframe = data.frame(fread('/home/patrick/code/pschupp/Singleton-analyses/mike-grant-endothelial-pseudobulk-analysis/SyntheticDatasets/SyntheticDataset1_25pcntCells_50pcntVar_100samples_08-25-51.csv'))
rnaDataframe$x = make.unique(rnaDataframe$x)
source('~/code/oldham-lab/FindModules/FindModules.lint.par.R')
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
wilcox.test(dat[,which(clusters$clust == 'Endothelial cells')+1], dat[,-(which(clusters$clust == 'Endothelial cells')+1)])
wilcox.test(dat[ind], x[!ind], alternative='greater')$p.value

```
