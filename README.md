# ChromothripsisPaperFigures
Scripts to used for the epigenetic data analysis for the chromothripsis paper.

## Log2 fold change comparing 117 and 152

## Normalised read counts plot comparing 117 and 152

## Distance plots 

## Nearest SV to differential genes

### For differential expression

We were interested in determining the nearest SV to each differential and non-differential gene so that it was possible to determine if there was a relationship between SV presence and differential genes. To determine this we used nearestSVisoseq.py:

python nearestSVisoseq.py

> python nearestSVisoseq.py  \
> --chromothriptic_bedpe data/chromothriptic_SVs_hg38_nogermline.bedpe \
> --wild_type_bedpe data/wild_type_SVs_hg38_nogermline.bedpe \
> --diff_genes data/differentialgenes.csv \
> --nondiff_genes data/nonDifferentialgenes.csv \
> --ensemble_info data/ensembleGeneID.csv
> --outputdir outs

This script requires pandas 

### For differential peaks from ChIP-seq and ATAC-seq Peaks

### plot distance to nearest SV

We wanted to see how far SV where from genes and whether this affected whether the genes were differential or not. For this we used plotNearestSVdistanceIsoseq.py:

python plotNearestSVdistanceIsoseq.py --inputdir outs --outputdir outs

You need seaborn, scipy, math, pandas and numpy.

### plot SV size to nearest SV

We wanted to see which SV types and sizes were nearest to the to the differential genes. For this we used plotNearestSVsizeDistributionIsoseq.py:

python plotNearestSVsizeDistributionIsoseq.py --inputdir outs/ --outputdir outs/

You need matplotlib, scipy, pandas and numpy.
