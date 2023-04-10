#from numpy.core.numeric import cross
from pydoc import doc
import site
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math 
import numpy as np 
plt.rcParams['pdf.fonttype'] = 42
from scipy import stats
import argparse 

parser = argparse.ArgumentParser(description='usage:  allelicFoldChange_Chip_ATAC_117_152.py   --outputdir outputdir')
parser.add_argument('--outputdir', help='path to output directory',required=True)
args = parser.parse_args()

outputdir=args.outputdir
sns.set()
sns.set_style(style='white') 
fig, axes = plt.subplots(1, 3, sharex=False, figsize=(9,10))

chromosomes=["chr1","chr6","chr9"]

color_dict = dict({'WTSI-OESO_117':'powderblue',
                  'WTSI-OESO_152':"lightgrey"})


for chr in range(len(chromosomes)):
    chromosome=chromosomes[chr]
    #set up the order of samples
    if chromosomes[chr]=="chr1":
        order=["WTSI-OESO_117","WTSI-OESO_152"]
        labels=["117","152"]
        correct=["no","no"]
        allele=[0,0]
    elif chromosomes[chr]=="chr6":
        order=["WTSI-OESO_117","WTSI-OESO_152"]
        labels=["117","152"]
        correct=["no","no"]
        allele=[0,0]
    elif chromosomes[chr]=="chr9":
        order=["WTSI-OESO_117","WTSI-OESO_152"]
        labels=["117","152"]
        correct=["yes","no"]
        allele=[2,0]
    input_files=[]
    for sample in order:
        input_files.append("{}{}{}{}{}".format("peaks_data/",sample,"_",chromosome,"_filtered_final_recorrected.csv"))
    #add the differential gene names to the dictionary
    peaks = pd.read_csv(input_files[0],sep="\t",header=None,skiprows=0)
    peaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr"]
    peak_ids=list(peaks["chr"].astype(str)+"_"+peaks["start"].astype(str)+"_"+peaks["end"].astype(str))
    peaks_dict={x:{} for x in peak_ids} #make the contig keys
    for sample in range(len(order)):
        foldchange = pd.read_csv(input_files[sample],sep="\t",header=None,skiprows=1)
        #if differential == "differential":
        foldchange.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr"]
        foldchange_ids=list(foldchange["chr"].astype(str)+"_"+foldchange["start"].astype(str)+"_"+foldchange["end"].astype(str))
        for item in range(len(foldchange_ids)):
            for x in peaks_dict:
                if foldchange_ids[item].split("_")[0]==x.split("_")[0] and int(foldchange_ids[item].split("_")[-2])-1000 < int(x.split("_")[-2]) < int(foldchange_ids[item].split("_")[-2])+1000 and int(foldchange_ids[item].split("_")[-1])-1000 < int(x.split("_")[-1]) < int(foldchange_ids[item].split("_")[-1])+1000:
                    #correct for copy number
                    if foldchange["new_fdr"][item]>0.05:
                        diff="nondiff"
                    else:
                        diff="diff"
                    fold="NA"
                    #for alleles which have while copy number gains, you need to half the binding values
                    if allele[sample]==0:
                        ct = 2**foldchange["Conc_chromothr"][item]
                        wt = 2**foldchange["Conc_wildtype"][item]
                    elif allele[sample]==1:
                        ct = (2**foldchange["Conc_chromothr"][item])/2
                        wt = 2**foldchange["Conc_wildtype"][item]
                    elif allele[sample]==2:
                        ct = 2**foldchange["Conc_chromothr"][item]
                        wt = (2**foldchange["Conc_wildtype"][item])/2
                    #calculate fold change - a positive fold change is up on CT
                    if ct/wt==0:
                        fold=0
                    else:
                        fold=math.log2(ct/wt)
                    peaks_dict[x][order[sample]]={"reported_fold_change":foldchange["Fold"][item], "fold_change":fold,"colour":diff}
    #get the log fold change 
    sample_fold_change=[]
    sample_colours=[]
    sample_position=[]
    for item in peaks_dict.keys():
        if  len(peaks_dict[item])==2: 
            sample_fold_change.append(peaks_dict[item]['WTSI-OESO_117']['fold_change'])
            sample_colours.append('WTSI-OESO_117')
            sample_position.append(0)
            sample_fold_change.append(peaks_dict[item]['WTSI-OESO_152']['fold_change'])
            sample_colours.append('WTSI-OESO_152')
            sample_position.append(1)
    #make a dataframe of the data to plot
    df = pd.DataFrame({"sample":sample_position, "fold_change":sample_fold_change, "colour":sample_colours})
    #append to plot
    sns.violinplot(ax=axes[chr], x="sample", y="fold_change", data=df, hue ="colour",palette=color_dict,scale="count")
    #sns.stripplot(ax=axes[chr], x="sample", y="fold_change", hue ="colour", data=df, palette=color_dict, jitter=True, s=4)
    axes[chr].set_title(chromosomes[chr])
    axes[chr].set_xlabel("") # same for y axis.
    #ind = np.arange(2)  # the x locations for the groups
    #axes[chr].set_xticks(ind, labels)
    axes[chr].set_xticklabels(labels)#, rotation=70)
    axes[chr].legend([],[], frameon=False)
    axes[chr].set_ylim(-15,17)
    axes[chr].axhline(y=0, color='grey', linestyle='-')
    #axes[x_pos,y_pos].set_xlim(-0.5,3.5)
    height=13.5
    i=0 
    for sample_number in range(len(order)):
        if sample_number!=len(order)-1:
            for number in range(sample_number+1,len(order)):
                i=i+1
                print(sample_number,number)
                statistic =  "p value = " + "{:.2e}".format(stats.wilcoxon(list(df.loc[df['sample'] == sample_number]["fold_change"]),list(df.loc[df['sample'] == number]["fold_change"]))[1], precision=2, unique=False, fractional=False, trim='k')
                axes[chr].plot([sample_number, sample_number, number, number], [height+i+i+0.1, height+i+i+0.35, height+i+i+0.35, height+i+i+0.1], lw=1.5, c='k')
                axes[chr].text((sample_number+number)*.5, height+i+i+0.35, statistic, ha='center', va='bottom',fontsize=13)
    #axes[x_pos,y_pos].set_xlim(-0.5,3.5)
    plt.tight_layout()



#axes[2].legend(loc=("upper right"))
plt.savefig("{}{}".format(outputdir,'/Diff_foldchange_117_152.pdf')  )

