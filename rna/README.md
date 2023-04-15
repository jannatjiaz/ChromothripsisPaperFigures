# Chromothripsis paper figure 5

This README provides scripts used to make figures used in figure 5. DESeq2 was used to proccess the haplotype resovled bam files containing reads which were explicity assigned to one haplotype over the other. After differentially expressed genes were called, those that were in regions of copy number variation or loss of hetrozygosity were removed as these were differential due to copy number effects. Genes which were fragmented by an SV were also removed. The remaining genes were explored further...

## Log2 fold change comparing 117 and 152

Filtered outputs from DESeq2 can be fed directly into allelicFoldChange_isoseq_117_152.py to see how the distribution of differences between alleles varies in two samples derived from the same donor. Use:

> python allelicFoldChange_isoseq_117_152.py --outputdir outs/

These scirpts require pandas, seaborn, matplotlib, numpy and scipy

## Normalised read counts plot comparing 117 and 152

DESeq2 can output the normalised count martrix and this can be used to determine how expression changes on the alleles over time. This direct alleic comparison can be done by:

> python directAllelicComparison_isoseq_117_152.py  --outputdir outs/

## Nearest SV to differential genes

We were interested in investigating the relationship between differential genes and differential peaks. For this peaks were split into active and inactive peaks and the distance to genes with a higher, lower and non-differential expression was dertermined using:

> python distanceofgenestopeaks.py
