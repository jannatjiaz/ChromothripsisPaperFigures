
import pandas as pd
import numpy as np
import argparse
import math
from collections import Counter
import matplotlib.pyplot as plt
from pylab import *
import scipy.stats as stats
import seaborn as sns

plt.rcParams['pdf.fonttype'] = 42


parser = argparse.ArgumentParser(description='usage:  plot_distanceofgenestopeaks.py  --differentialPeaksInput differentialPeaksInput --nonDifferentialPeaksInput nonDifferentialPeaksInput --activity active --on_CT on_CT')
parser.add_argument('--differentialPeaksInput', help='path to differential peaks file on the chromosome with peak',required=True)
parser.add_argument('--nonDifferentialPeaksInput', help='path to non-differential peaks file  on the chromosome with peak',required=True)
parser.add_argument('--activity', help='active or inactive depending of which peaks you are looking at',required=True)
parser.add_argument('--on_CT', help='higher or lower depending on whether you want to look at peaks that are stronger on the chromothriptic or weaker on the chromothriptic ',required=True)
args = parser.parse_args()

#start with genes that are higher on CT

differentialPeaksInput=args.differentialPeaksInput
nonDifferentialPeaksInput=args.nonDifferentialPeaksInput
activity=args.activity
on_CT=args.on_CT



#start with genes that are higher on CT

differentialPeaks = pd.read_csv(differentialPeaksdifferentialPeaks,sep="\t",header=None, skiprows=0)
differentialPeaks.columns = ["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr","distance_to_diffpeak_middle","distance_to_nondiffpeak_middle"]

nonDifferentialPeaks = pd.read_csv(nonDifferentialPeaksInput,sep="\t",header=None, skiprows=0)
nonDifferentialPeaks.columns = ["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr","distance_to_diffpeak_middle","distance_to_nondiffpeak_middle"]


#binned differnce plot for differential peaks
#for differential genes
diff_distance_differential = []
nondiff_distance_differential = []
if on_CT == "higher":
    for item in range(len(differentialPeaks)):
        #get the absolute value of the difference
        if differentialPeaks["log2FoldChange"][item]<0:
            distance = np.log10(differentialPeaks["distance_to_diffpeak_middle"][item])
            diff_distance_differential.append(distance)
            distance = np.log10(differentialPeaks["distance_to_nondiffpeak_middle"][item])
            nondiff_distance_differential.append(distance)
else:
    for item in range(len(differentialPeaks)):
        #get the absolute value of the difference
        if differentialPeaks["log2FoldChange"][item]>0:
            distance = np.log10(differentialPeaks["distance_to_diffpeak_middle"][item])
            diff_distance_differential.append(distance)
            distance = np.log10(differentialPeaks["distance_to_nondiffpeak_middle"][item])
            nondiff_distance_differential.append(distance)

#for nondifferential genes
diff_distance_nondifferential = []
nondiff_distance_nondifferential = []
if on_CT == "higher":
    for item in range(len(nonDifferentialPeaks)):
        #get the absolute value of the difference
        distance = np.log10(nonDifferentialPeaks["distance_to_diffpeak_middle"][item])
        diff_distance_nondifferential.append(distance)
        distance = np.log10(nonDifferentialPeaks["distance_to_nondiffpeak_middle"][item])
        nondiff_distance_nondifferential.append(distance)
else:
    for item in range(len(nonDifferentialPeaks)):
        #get the absolute value of the difference
        distance = np.log10(nonDifferentialPeaks["distance_to_diffpeak_middle"][item])
        diff_distance_nondifferential.append(distance)
        distance = np.log10(nonDifferentialPeaks["distance_to_nondiffpeak_middle"][item])
        nondiff_distance_nondifferential.append(distance)


#turn it into an array
all_distances = []
#chromothriptic
diff_distance_differential = np.array(diff_distance_differential)
all_distances.append(diff_distance_differential)
nondiff_distance_differential = np.array(nondiff_distance_differential)
all_distances.append(nondiff_distance_differential)
#wildtype
diff_distance_nondifferential = np.array(diff_distance_nondifferential)
all_distances.append(diff_distance_nondifferential)
nondiff_distance_nondifferential = np.array(nondiff_distance_nondifferential)
all_distances.append(nondiff_distance_nondifferential)

#plot a scatter plot overlayed with a boxplot

plt.figure(figsize =(10, 7))
sns.set(style='ticks', context='talk')

if on_CT=="higher":
    sns.distplot(all_distances[1],color='skyblue',hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  label = "distance to non-differential peaks")
    sns.distplot(all_distances[0], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='firebrick', label = "distance to differential peaks")
    plt.legend(labels=['distance to non-differential peaks','distance to differential peaks'])
else:
    sns.distplot(all_distances[1], color='skyblue',hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  label = "distance to non-differential peaks")
    sns.distplot(all_distances[0], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},   color='firebrick', label = "distance to differential peaks")
    plt.legend(labels=['distance to non-differential peaks','distance to differential peaks'])

#plt.ylim(0,3)
plt.xlim(1,8)

plt.xlabel('Distance (log10(bp))')
if on_CT=="higher":
    plt.title('Higher gene expression on chromothriptic allele')
else:
    plt.title('Higher gene expression on wild type allele')

non_diff_genes_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(nondiff_distance_nondifferential, diff_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')
diff_genes_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(diff_distance_differential, nondiff_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')


plt.figtext(0.25, 0.8, diff_genes_stats, ha="center",  wrap=True, horizontalalignment='center', fontsize=12)


plt.savefig("{}{}{}".format('nearestSVlogscale',on_CT,'_filtered_differential_genes_activepeaks.pdf')  )

plt.show()




plt.figure(figsize =(10, 7))
sns.set(style='ticks', context='talk')

if on_CT=="higher":
    sns.distplot(all_distances[3], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},color='skyblue', label = "distance to non-differential peaks")
    sns.distplot(all_distances[2], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3},  color='firebrick', label = "distance to differential peaks")
    plt.legend(labels=['distance to non-differential peaks','distance to differential peaks'])
else:
    sns.distplot(all_distances[3], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='skyblue', label = "distance to non-differential peaks")
    sns.distplot(all_distances[2], hist = False, kde = True, kde_kws = {'shade': True, 'linewidth': 3}, color='firebrick', label = "distance to differential peaks")
    plt.legend(labels=['distance to non-differential peaks','distance to differential peaks'])


#plt.ylim(0,3)
plt.xlim(1,8)

plt.xlabel('Distance (log10(bp))')
if on_CT=="higher":
    plt.title('Higher gene expression on chromothriptic allele')
else:
    plt.title('Higher gene expression on wild type allele')

non_diff_genes_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(nondiff_distance_nondifferential, diff_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')
diff_genes_stats = "p value = "+np.format_float_positional(stats.mannwhitneyu(diff_distance_differential, nondiff_distance_nondifferential, alternative='two-sided')[1], precision=2, unique=False, fractional=False, trim='k')

plt.figtext(0.25, 0.8, non_diff_genes_stats, ha="center",  wrap=True, horizontalalignment='center', fontsize=12)


plt.savefig("{}{}{}".format('nearestSVlogscale',on_CT,'_filtered_non_differential_genes_activepeaks.pdf')  )

plt.show()

