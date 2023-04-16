# Chromothripsis paper figure 5

This README provides scripts used to make figures used in figure 5. DESeq2 was used to proccess the haplotype resovled bam files containing reads which were explicity assigned to one haplotype over the other. After differentially expressed genes were called, those that were in regions of copy number variation or loss of hetrozygosity were removed as these were differential due to copy number effects. Genes which were fragmented by an SV were also removed. The remaining genes were explored further...

## Log2 fold change comparing 117 and 152

Filtered outputs from DESeq2 can be fed directly into allelicFoldChange_isoseq_117_152.py to see how the distribution of differences between alleles varies in two samples derived from the same donor. Use:

> python allelicFoldChange_isoseq_117_152.py --outputdir outs/

This scirpt require pandas, seaborn, matplotlib, numpy and scipy

## Normalised read counts plot comparing 117 and 152

DESeq2 can output the normalised count martrix and this can be used to determine how expression changes on the alleles over time. This direct alleic comparison can be done by:

> python directAllelicComparison_isoseq_117_152.py  --outputdir outs/

## Nearest peak to genes analysis

We wanted to investigate the realationship between chromatin accessibility peaks and genes. For this we determined distance of differential and non-differentail genes to peaks. For example, for differential genes:

> python distanceofgenestopeaks.py  --differential differential \
> --atacchip_diff example_inputs/nearestSV_differentialPeaks_noCNchange_noSVinteracting.csv \
> --atacchip_nondiff example_inputs/nearestSV_nonDifferentialPeaks_noCNchange_noSVinteracting.csv \
> --genes_input_dir example_inputs --genes_file_name differentialgenes_noCNchange_noSVinteracting.csv 
> --ensemble_input example_inputs/ensembleGeneID.csv

This scirpt require pandas.

From the output of diffbind you will need to associate peaks with active or repressive marks or CTCF. This can be done using bedtools, for example:

> bedtools closest -a nearestSV_nonDifferentialPeaks.csv -b <(cat macs2/chromothriptic_reads_*_assigned_peaks.narrowPeak macs2/chromothriptic_reads_*_assigned_peaks.broadPeak | cut -f 1,2,3,4) | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,17 | sort | uniq > nearestSV_nonDifferentialPeaks_to_inactive_histone_mark.csv



Once we have determined the distane between peaks and gene, the peaks were split into active and inactive peaks and the distance to genes with a higher, lower and non-differential expression was dertermined using:

> python plotNearestSv_chipatac_distplot.py  --differentialPeaksInput nearestSV_differentialPeaks_to_active_histone_mark.csv \
> --nonDifferentialPeaksInput nearestSV_nonDifferentialPeaks_to_active_histone_mark.csv --activity active --on_CT higher

This script requires pandas, numpy, math, collections, matplotlib, pylab, scipy and seaborn. 

