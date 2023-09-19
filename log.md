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
[ ] Perform pseudobulk
[ ] Do FM of pseudobulk
[ ] Perform enrichments of networks
[ ] Perform DE in snRNA-seq

## Log

Begin documenting work.

## Pseudobulk

```{.r}

