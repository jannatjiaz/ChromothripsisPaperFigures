# Chromothripsis paper figures 4

This README provides scripts used to make figures used in figure 4. Diffbind was used to proccess the haplotype resovled bam files containing reads which were explicity assigned to one haplotype over the other. After differential peaks were called, those that were in regions of copy number variation or loss of hetrozygosity were removed as these were differential due to copy number effects. Peaks which were directly bisected by an SV were also removed. The remaining peaks were explored further... 

## Log2 fold change comparing 117 and 152

Filtered outputs from DiffBind can be fed directly into allelicFoldChange_Chip_ATAC_117_152.py to see how the distribution of differences between alleles varies in two samples derived from the same donor. Use:

> python allelicFoldChange_Chip_ATAC_117_152.py --outputdir outs/

This scirpt require pandas, seaborn, matplotlib, numpy and scipy

## Normalised read counts plot comparing 117 and 152

In order to compare the ChIP and ATAC-seq peaks, samtools was used to determine the coverage at each peak called in the ChIP and ATAC data. For each peak, the number of reads that were present in that region was determined for all ChIP-seq marks and for the ATAC-seq runs. This produced a file where each row was a peak and the columns were representing the counts of each mark at the peak site. From this a direct allelic comparison can be done:

> python directAllelicComparison_Chip_ATAC_117_152.py  --outputdir outs/


## Nearest SV to differential genes and peaks

We were interested in determining the nearest SV to each differential and non-differential gene and peak so that it was possible to determine if there was a relationship between SV presence and differential peaks. 


To determine this for ChIP-seq and ATAC-seq peaks, we used nearestSVpeaks.py:
> python nearestSVpeaks.py  \
> --chromothriptic_bedpe example_inputs/chromothriptic_SVs_hg38_nogermline.bedpe \
> --wild_type_bedpe example_inputs/wild_type_SVs_hg38_nogermline.bedpe \
> --diff_peaks example_inputs/differentialPeaks.csv \
> --nondiff_peaks example_inputs/nonDifferentialPeaks.csv \
> --outputdir outs

This script requires pandas. The output of this was filtered to remove regions with copy number altereations.

### plot distance to nearest SV

We wanted to visulaise how far SV where from genes and whether this affected whether the genes were differential or not. For this we used plotNearestSVdistancePeaks.py. For example:

> python plotNearestSv_chipatac_distplot.py  --differentialPeaksInput example_inputs/nearestSV_differentialPeaks_noCNchange_noSVinteracting.csv \
> --nonDifferentialPeaksInput example_inputs/nearestSV_nonDifferentialPeaks_noCNchange_noSVinteracting.csv \
> --differentialPeaksOtherInput example_inputs/nearestSV_differential_on_other_allele_Peaks_noCNchange_noSVinteracting.csv \
> --nonDifferentialPeaksOtherInput example_inputs/nearestSV_nonDifferential_on_other_allele_Peaks_noCNchange_noSVinteracting.csv --on_CT higher

You need seaborn math, pandas, numpy, collections, matplotlib, pylab, scipy and seaborn.

### plot SV size of nearest SV (and type)

We wanted to see which what the effect of distance to SV had on the fold change difference in peak height. You need a peak file with a column that annotates which peak is assocated which which histone modification. For this we used plotNearestSv_chipatac_distplot.py:

> python distance_foldchange_effect.py  --differentialPeaksInput example_inputs/nearestSV_differentialPeaks_noCNchange_noSVinteracting_withpeakinfo.csv \
> --nonDifferentialPeaksInput example_inputs/nearestSV_nonDifferentialPeaks_noCNchange_noSVinteracting_withpeakinfo.csv \
> --differentialPeaksOtherInput example_inputs/nearestSV_differential_on_other_allele_Peaks_noCNchange_noSVinteracting_withpeakinfo.csv \
> --nonDifferentialPeaksOtherInput example_inputs/nearestSV_nonDifferential_on_other_allele_Peaks_noCNchange_noSVinteracting_withpeakinfo.csv

This scirpt needs matplotlib, math, pandas and numpy.
