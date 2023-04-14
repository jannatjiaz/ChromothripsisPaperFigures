
import pandas as pd
import numpy as np
import math
import argparse
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

### CHIP/ATAC

parser = argparse.ArgumentParser(description='usage:  plotNearestSv_chipatac_distplot.py  --differentialPeaksInput differentialPeaksInput --nonDifferentialPeaksInput nonDifferentialPeaksInput \
    --differentialPeaksOtherInput differentialPeaksOtherInput --nonDifferentialPeaksOtherInput nonDifferentialPeaksOtherInput')
parser.add_argument('--differentialPeaksInput', help='path to differential peaks file on the chromosome with peak',required=True)
parser.add_argument('--nonDifferentialPeaksInput', help='path to non-differential peaks file  on the chromosome with peak',required=True)
parser.add_argument('--differentialPeaksOtherInput', help='path to differential peaks file on the chromosome without peak',required=True)
parser.add_argument('--nonDifferentialPeaksOtherInput', help='path to non-differential peaks file on the chromosome without peak',required=True)
args = parser.parse_args()


differentialPeaksInput = "{}".format("/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/hifiasm_assembly/isoseq_hifiasm_assembly/DESeq2/new/tmp/actachip/split_by_mark_final/nearestSV_differentialPeaks_noCNchange_noSVinteracting_peak.csv")
nonDifferentialPeaksInput = "{}".format("/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/hifiasm_assembly/isoseq_hifiasm_assembly/DESeq2/new/tmp/actachip/split_by_mark_final/nearestSV_nonDifferentialPeaks_noCNchange_noSVinteracting_peak.csv")
differentialPeaksOtherInput = "{}".format("/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/hifiasm_assembly/isoseq_hifiasm_assembly/DESeq2/new/tmp/actachip/split_by_mark_final/nearestSV_differential_on_other_allele_Peaks_noCNchange_noSVinteracting_peak.csv")
nonDifferentialPeaksOtherInput = "{}".format("/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/hifiasm_assembly/isoseq_hifiasm_assembly/DESeq2/new/tmp/actachip/split_by_mark_final/nearestSV_nonDifferential_on_other_allele_Peaks_noCNchange_noSVinteracting_peak.csv")


color_dict={}
color_dict["CTCF"]="black"
color_dict["H3K27ac"]="orange"
color_dict["H3K27me3"]="crimson"
color_dict["H3K4me3"]="dodgerblue"

differentialPeaks = pd.read_csv(differentialPeaksInput,sep="\t",header=None, skiprows=1)
differentialPeaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle","peak"]
nonDifferentialPeaks = pd.read_csv(nonDifferentialPeaksInput,sep="\t",header=None, skiprows=1)
nonDifferentialPeaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct","peaks"]
differentialPeaksOther = pd.read_csv(differentialPeaksOtherInput,sep="\t",header=None, skiprows=1)
differentialPeaksOther.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle","peak"]
nonDifferentialPeaksOther = pd.read_csv(nonDifferentialPeaksOtherInput,sep="\t",header=None, skiprows=1)
nonDifferentialPeaksOther.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct","peak"]


distance_stronger=[]
fold_change_stronger=[]
colour_stronger=[]

distance_weaker=[]
fold_change_weaker=[]
colour_weaker=[]


for item in range(len(differentialPeaks)):
    #get the absolute value of the difference
    if differentialPeaks["Fold"][item]>0:
        distance_val = np.log10(differentialPeaks["SV_distance_to_peak_middle"][item])
        distance_stronger.append(distance_val)
        colour_stronger.append(color_dict[differentialPeaks['peak'][item].split("_")[5] ])
        fold_change_stronger.append(differentialPeaks["Fold"][item])
    else:
        distance_val = np.log10(differentialPeaksOther["SV_distance_to_peak_middle"][item])
        distance_weaker.append(distance_val)
        colour_weaker.append(color_dict[differentialPeaks['peak'][item].split("_")[5]])
        fold_change_weaker.append(-1*differentialPeaksOther["Fold"][item])

for item in range(len(nonDifferentialPeaks)):
    if nonDifferentialPeaks["Fold"][item]>0:
        distance_stronger.append(math.log10(nonDifferentialPeaks["SV_distance_to_peak_middle_ct"][item]))
        colour_stronger.append("lightgrey")
        fold_change_stronger.append(nonDifferentialPeaks["Fold"][item])
    else:
        distance_weaker.append(math.log10(nonDifferentialPeaks["SV_distance_to_peak_middle_ct"][item]))
        colour_weaker.append("lightgrey")
        fold_change_weaker.append(-1*nonDifferentialPeaks["Fold"][item])


plt.scatter(distance_stronger,fold_change_stronger, color=colour_stronger, s=7, alpha=1)
plt.title("Stonger")
plt.xlabel("distance")
plt.ylabel("log2fold change in peak height expression")
plt.xlim(1,7)
plt.ylim(-0.5,14)
plt.savefig("{}".format('allPeaks_distance_foldchange_strongerbinding.pdf')  )

plt.show()


plt.scatter(distance_weaker,fold_change_weaker, color=colour_weaker, s=7, alpha=1)
plt.title("Weaker")
plt.xlabel("distance")
plt.ylabel("log2fold change in peak height expression")
plt.xlim(1,7)
plt.ylim(-0.5,14)
plt.savefig("{}".format('allPeaks_distance_foldchange_weakerbinding.pdf')  )

plt.show()
