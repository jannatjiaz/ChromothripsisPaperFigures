# Chromothripsis paper figure 5

This README provides scripts used to make figures used in figure 5. DESeq2 was used to proccess the aligned haplotype resovled bam files and the outputs from these were used for plotting puposes. 

## Log2 fold change comparing 117 and 152

After filtering the peaks to remove copy number altered regions, outputs from DESeq2 can be fed directly into allelicFoldChange_isoseq_117_152.py to see how the distribution of differences between alleles varies in two samples derived from the same donor. Use:

> python allelicFoldChange_isoseq_117_152.py --outputdir outs/

These scirpts require pandas, seaborn, matplotlib, numpy and scipy

## Normalised read counts plot comparing 117 and 152

DESeq2 can output the normalised count martrix and this can be used to determine how expression changes on the alleles over time. This direct alleic comparison can be done by:

> python directAllelicComparison_isoseq_117_152.py  --outputdir outs/

## Nearest SV to differential genes

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

### plot distance to nearest SV

We wanted to visulaise how far SV where from genes and whether this affected whether the genes were differential or not. For this we used plotNearestSVdistanceIsoseq.py

> python plotNearestSVdistanceIsoseq.py --inputdir outs --outputdir outs

You need seaborn, scipy, math, pandas and numpy.

### plot SV size of nearest SV (and type)

We wanted to see which SV types and sizes were nearest to the to the differential genes and peaks. For this we used plotNearestSVsizeDistributionIsoseq.py:

> python plotNearestSVsizeDistributionIsoseq.py --inputdir outs/ --outputdir outs/

You need matplotlib, scipy, pandas and numpy.
