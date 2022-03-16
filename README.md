# ChromothripsisPaperFigures
Scripts to used for the epigenetic data analysis for the chromothripsis paper.

## Log2 fold change comparing 117 and 152

## Normalised read counts plot comparing 117 and 152

## Nearest SV to differential genes and peaks

We were interested in determining the nearest SV to each differential and non-differential gene and peak so that it was possible to determine if there was a relationship between SV presence and differential genes or peaks. 

To determine this for genes, we used nearestSVisoseq.py:

> python nearestSVisoseq.py  \
> --chromothriptic_bedpe data/chromothriptic_SVs_hg38_nogermline.bedpe \
> --wild_type_bedpe data/wild_type_SVs_hg38_nogermline.bedpe \
> --diff_genes data/differentialgenes.csv \
> --nondiff_genes data/nonDifferentialgenes.csv \
> --ensemble_info data/ensembleGeneID.csv
> --outputdir outs

This script requires pandas 

To determine this for ChIP-seq and ATAC-seq peaks, we used nearestSVpeaks.py:
> python nearestSVpeaks.py  \
> --chromothriptic_bedpe data/chromothriptic_SVs_hg38_nogermline.bedpe \
> --wild_type_bedpe data/wild_type_SVs_hg38_nogermline.bedpe \
> --diff_peaks data/differentialPeaks.csv \
> --nondiff_peaks data/nonDifferentialPeaks.csv \
> --outputdir outs

This script requires pandas 

### plot distance to nearest SV

We wanted to visulaise how far SV where from genes and whether this affected whether the genes were differential or not. For this we used plotNearestSVdistanceIsoseq.py and plotNearestSVdistancePeaks.py:

python plotNearestSVdistanceIsoseq.py --inputdir outs --outputdir outs
python plotNearestSVdistancePeaks.py --inputdir outs --outputdir outs

You need seaborn, scipy, math, pandas and numpy.

### plot SV size of nearest SV (and type)

We wanted to see which SV types and sizes were nearest to the to the differential genes and peaks. For this we used plotNearestSVsizeDistributionIsoseq.py and plotNearestSVsizeDistributionPekas.py:

python plotNearestSVsizeDistributionIsoseq.py --inputdir outs/ --outputdir outs/
python plotNearestSVsizeDistributionPekas.py --inputdir outs/ --outputdir outs/

You need matplotlib, scipy, pandas and numpy.
