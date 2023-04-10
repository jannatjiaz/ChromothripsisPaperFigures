# Chromothripsis paper figures 4

This README provides scripts used to make figures used in figure 4. Diffbind was used to proccess the aligned haplotype resovled bam files and the outputs from these were used for plotting puposes. 

## Log2 fold change comparing 117 and 152

After filtering the peaks to remove copy number altered regions, outputs from DiffBind can be fed directly into allelicFoldChange_Chip_ATAC_117_152.py to see how the distribution of differences between alleles varies in two samples derived from the same donor. Use:

> python allelicFoldChange_Chip_ATAC_117_152.py --outputdir outs/

These scirpts require pandas, seaborn, matplotlib, numpy and scipy

## Normalised read counts plot comparing 117 and 152

In order to compare the ChIP and ATAC-seq peaks, samtools was used to determine the coverage at each peak called in the ChIP and ATAC data. For each peak, the number of reads that were present in that region was determined for all ChIP-seq marks and for the ATAC-seq runs. This produced a file where each row was a peak and the columns were representing the counts of each mark at the peak site. From this a direct allelic comparison can be done:

> python directAllelicComparison_Chip_ATAC_117_152.py  --outputdir outs/


## Nearest SV to differential genes and peaks

We were interested in determining the nearest SV to each differential and non-differential gene and peak so that it was possible to determine if there was a relationship between SV presence and differential peaks. 


To determine this for ChIP-seq and ATAC-seq peaks, we used nearestSVpeaks.py:
> python nearestSVpeaks.py  \
> --chromothriptic_bedpe data/chromothriptic_SVs_hg38_nogermline.bedpe \
> --wild_type_bedpe data/wild_type_SVs_hg38_nogermline.bedpe \
> --diff_peaks data/differentialPeaks.csv \
> --nondiff_peaks data/nonDifferentialPeaks.csv \
> --outputdir outs

This script requires pandas 

### plot distance to nearest SV

We wanted to visulaise how far SV where from genes and whether this affected whether the genes were differential or not. For this we used plotNearestSVdistancePeaks.py:

> python plotNearestSVdistancePeaks.py --inputdir outs --outputdir outs

You need seaborn, scipy, math, pandas and numpy.

### plot SV size of nearest SV (and type)

We wanted to see which SV types and sizes were nearest to the to the differential genes and peaks. For this we used plotNearestSVsizeDistributionPekas.py:

> python plotNearestSVsizeDistributionPekas.py --inputdir outs/ --outputdir outs/

You need matplotlib, scipy, pandas and numpy.
