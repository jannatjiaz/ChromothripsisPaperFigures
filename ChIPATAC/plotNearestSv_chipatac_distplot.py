import pandas as pd
import numpy as np 
import argparse 
import math
from collections import Counter
import matplotlib.pyplot as plt
from pylab import *
import scipy.stats as stats
import seaborn as sns
import argparse

plt.rcParams['pdf.fonttype'] = 42


parser = argparse.ArgumentParser(description='usage:  plotnearestSVpeaks.py  --chromothriptic_bedpe chromothriptic_bedpe --wild_type_bedpe wild_type_bedpe \
    --diff_genes diff_genes_input -nondiff_genes nondiff_genes_input --ensemble_info ensemble_input')
parser.add_argument('--differentialPeaksInput', help='path to differential peaks file on the chromosome with peak',required=True)
parser.add_argument('--nonDifferentialPeaksInput', help='path to non-differential peaks file  on the chromosome with peak',required=True)
parser.add_argument('--differentialPeaksOtherInput', help='path to differential peaks file on the chromosome with peak',required=True)
parser.add_argument('--nonDifferentialPeaksOtherInput', help='path to non-differential peaks file  on the chromosome with peak',required=True)
parser.add_argument('--on_CT', help='higher or lower depending on whether you want to look at peaks that are stronger on the chromothriptic or weaker on the chromothriptic ',required=True)
args = parser.parse_args()

#start with genes that are higher on CT 

differentialPeaksInput=args.differentialPeaksInput
nonDifferentialPeaksInput=args.nonDifferentialPeaksInput
differentialPeaksOtherInput=args.differentialPeaksOtherInput
nonDifferentialPeaksOtherInput=args.nonDifferentialPeaksOtherInput
on_CT=args.on_CT

differentialPeaks = pd.read_csv(differentialPeaksInput,sep="\t",header=None, skiprows=1)
differentialPeaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle"]
nonDifferentialPeaks = pd.read_csv(nonDifferentialPeaksInput,sep="\t",header=None, skiprows=1)
nonDifferentialPeaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct"]

differentialPeaksOther = pd.read_csv(differentialPeaksOtherInput,sep="\t",header=None, skiprows=1)
differentialPeaksOther.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle"]
nonDifferentialPeaksOther = pd.read_csv(nonDifferentialPeaksOtherInput,sep="\t",header=None, skiprows=1)
nonDifferentialPeaksOther.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct"]

#binned differnce plot for differential peaks
wt_distance_differential = []
ct_distance_differential = []
if on_CT == "higher":
    for item in range(len(differentialPeaks)):
        #get the absolute value of the difference 
        if differentialPeaks["Conc_wildtype"][item]<differentialPeaks["Conc_chromothr"][item]:
            distance = np.log10(differentialPeaksOther["SV_distance_to_peak_middle"][item])
            wt_distance_differential.append(distance)
            distance = np.log10(differentialPeaks["SV_distance_to_peak_middle"][item])
            ct_distance_differential.append(distance)
else:
    for item in range(len(differentialPeaks)):
        #get the absolute value of the difference 
        if differentialPeaks["Conc_wildtype"][item]>differentialPeaks["Conc_chromothr"][item]:
            distance = np.log10(differentialPeaks["SV_distance_to_peak_middle"][item])
            wt_distance_differential.append(distance)
            distance = np.log10(differentialPeaksOther["SV_distance_to_peak_middle"][item])
            ct_distance_differential.append(distance)



wt_distance_nondifferential = []
ct_distance_nondifferential = []
if on_CT == "higher":
    for item in range(len(nonDifferentialPeaks)):
        #get the absolute value of the difference 
        #distance = math.ceil(nonDifferentialPeaks["SV_distance_to_peak_middle"][item] / 10000) * 10000
        if nonDifferentialPeaks["SV_distance_to_peak_middle_wt"][item]!=0:
            distance = np.log10(nonDifferentialPeaks["SV_distance_to_peak_middle_wt"][item])
            wt_distance_nondifferential.append(distance)
        else:
            wt_distance_nondifferential.append(0)
        if nonDifferentialPeaks["SV_distance_to_peak_middle_ct"][item]!=0:
            distance = np.log10(nonDifferentialPeaks["SV_distance_to_peak_middle_ct"][item])
            ct_distance_nondifferential.append(distance)
        else:
            ct_distance_nondifferential.append(0)
else:
    for item in range(len(nonDifferentialPeaks)):
        #get the absolute value of the difference 
        #distance = math.ceil(nonDifferentialPeaks["SV_distance_to_peak_middle"][item] / 10000) * 10000
        if nonDifferentialPeaksOther["SV_distance_to_peak_middle_wt"][item]!=0:
            distance = np.log10(nonDifferentialPeaksOther["SV_distance_to_peak_middle_ct"][item])
            wt_distance_nondifferential.append(distance)
        else:
            wt_distance_nondifferential.append(0)
        if nonDifferentialPeaksOther["SV_distance_to_peak_middle_ct"][item]!=0:
            distance = np.log10(nonDifferentialPeaksOther["SV_distance_to_peak_middle_wt"][item])
            ct_distance_nondifferential.append(distance)
        else:
            ct_distance_nondifferential.append(0)


#turn it into an array 
all_distances = []
#chromothriptic 
ct_distance_differential = np.array(ct_distance_differential)
all_distances.append(ct_distance_differential)
ct_distance_nondifferential = np.array(ct_distance_nondifferential)
all_distances.append(ct_distance_nondifferential)
#wildtype
wt_distance_differential = np.array(wt_distance_differential)
all_distances.append(wt_distance_differential)
wt_distance_nondifferential = np.array(wt_distance_nondifferential)
all_distances.append(wt_distance_nondifferential)

#plot a scatter plot overlayed with a boxplot

plt.figure(figsize =(10, 7))
sns.set(style='ticks', context='talk')

if on_CT=="higher":
    sns.distplot(all_distances[1],color='skyblue',hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  label = "equal on CT and WT")
    sns.distplot(all_distances[0], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='firebrick', label = "stronger on CT/Weaker on WT")
    plt.legend(labels=['equal on CT and WT','weaker on CT/Weaker on WT'])
else:
    sns.distplot(all_distances[1], color='skyblue',hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  label = "equal on CT and WT")
    sns.distplot(all_distances[0], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},   color='firebrick', label = "weaker on CT/stronger on WT")
    plt.legend(labels=['equal on CT and WT','weaker on CT/Weaker on WT'])

plt.ylim(0,1)

plt.xlabel('Distance (log10(bp))')
if on_CT=="higher":
    plt.title('Stonger ATAC-seq and ChIP-seq peaks on chromothriptic allele CT')
else:
    plt.title('Weaker ATAC-seq and ChIP-seq peaks on chromothriptic allele CT')

wt_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(wt_distance_differential, wt_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')
ct_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(ct_distance_differential, ct_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')


plt.figtext(0.25, 0.8, ct_stats, ha="center",  wrap=True, horizontalalignment='center', fontsize=12)


plt.savefig("{}{}{}".format('nearestSVlogscale',on_CT,'_filtered_CT.pdf')  )

plt.show()




plt.figure(figsize =(10, 7))
sns.set(style='ticks', context='talk')

if on_CT=="higher":
    sns.distplot(all_distances[3], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},color='skyblue', label = "equal on CT and WT")
    sns.distplot(all_distances[2], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  color='firebrick', label = "stronger on CT/Weaker on WT")
    plt.legend(labels=['equal on CT and WT','weaker on CT/Weaker on WT'])
else:
    sns.distplot(all_distances[3], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='skyblue', label = "equal on CT and WT")
    sns.distplot(all_distances[2], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='firebrick', label = "weaker on CT/stronger on WT")
    plt.legend(labels=['equal on CT and WT','weaker on CT/Weaker on WT'])

plt.ylim(0,1)

plt.xlabel('Distance (log10(bp))')
if on_CT=="higher":
    plt.title('Stonger ATAC-seq and ChIP-seq peaks on chromothriptic allele WT')
else:
    plt.title('Weaker ATAC-seq and ChIP-seq peaks on chromothriptic allele WT')

wt_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(wt_distance_differential, wt_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')
ct_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(ct_distance_differential, ct_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')


plt.figtext(0.25, 0.8, wt_stats, ha="center",  wrap=True, horizontalalignment='center', fontsize=12)


plt.savefig("{}{}{}".format('nearestSVlogscale',on_CT,'_filtered_WT.pdf')  )

plt.show()



#####
