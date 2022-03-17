
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math 
import numpy as np
from scipy import stats
plt.rcParams['pdf.fonttype'] = 42    
import argparse 

parser = argparse.ArgumentParser(description='usage:  allelicFoldChange_isoseq_117_152.py   --outputdir outputdir')
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
    elif chromosomes[chr]=="chr6":
        order=["WTSI-OESO_117","WTSI-OESO_152"]
        labels=["117","152"]
    elif chromosomes[chr]=="chr9":
        order=["WTSI-OESO_117","WTSI-OESO_152"]
        labels=["117","152"]
    #get input files paths
    input_files=[]
    for sample in order:
        s=sample.split("_")[1]
        input_files.append("{}{}{}{}{}".format("rna_data/differentialExpression_",chromosome,"_",s,"_filtered_final_recorrected.csv"))        
    #add the differential gene names to the dictionary
    genes = pd.read_csv(input_files[0],sep="\t",header=None,skiprows=1)
    genes.columns = ["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr"]
    genes_dict={x:{} for x in genes["gene"]} #make the contig keys
    #fill the dictionary in with the log2foldchange info
    for sample in range(len(order)):
        genes = pd.read_csv(input_files[sample],sep="\t",header=None,skiprows=1)
        genes.columns = ["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr"]
        for item in range(len(genes)):
            if genes["gene"][item] in genes_dict:
                if genes["new_fdr"][item]>0.05:
                    diff="nondiff"
                else:
                    diff="diff"
                fold = -1*genes["log2FoldChange"][item]                
                genes_dict[genes["gene"][item]][order[sample]]={'fold':fold, "colour":diff}
    #get the log fold change 
    sample_fold_change=[]
    sample_colours=[]
    sample_position=[]
    for item in genes_dict.keys():
        if  len(genes_dict[item])==2: 
            sample_fold_change.append(genes_dict[item]['WTSI-OESO_117']['fold'])
            #sample_colours.append(genes_dict[item][order[sample_number]]['colour'])
            sample_colours.append('WTSI-OESO_117')
            sample_position.append(0)
            sample_fold_change.append(genes_dict[item]['WTSI-OESO_152']['fold'])
            #sample_colours.append(genes_dict[item][order[sample_number]]['colour'])
            sample_colours.append('WTSI-OESO_152')
            sample_position.append(1)
    #make a dataframe of the data to plot
    df = pd.DataFrame({"sample":sample_position, "fold_change":sample_fold_change, "colour":sample_colours})
    #make the position paraments for the plot
    if chr==0 or chr==1 or chr==2:
        x_pos=0
    else:
        x_pos=1
    if chr==0 or chr==3:
        y_pos=0
    elif chr==1 or chr==4:
        y_pos=1    
    else:
        y_pos=2
    #append to plot
    sns.violinplot(ax=axes[chr], x="sample", y="fold_change", data=df, hue ="colour",palette=color_dict,scale="count")
    #sns.stripplot(ax=axes[chr], x="sample", y="fold_change", hue ="colour", data=df, palette=color_dict, jitter=True, s=4)
    axes[chr].set_title(chromosomes[chr])
    axes[chr].set_xlabel("") # same for y axis.
    #ind = np.arange(4)  # the x locations for the groups
    #axes[chr].set_xticks(ind, labels)
    axes[chr].set_xticklabels(labels)#, rotation=20)
    axes[chr].legend([],[], frameon=False)
    axes[chr].set_ylim(-15,12)
    axes[chr].axhline(y=0, color='grey', linestyle='-')
    height=8.5
    i=0 
    for sample_number in range(len(order)):
        if sample_number!=len(order)-1:
            for number in range(sample_number+1,len(order)):
                i=i+1
                print(sample_number,number)
                statistic =  "p value = " + "{:.2e}".format(stats.wilcoxon(list(df.loc[df['sample'] == 0]["fold_change"]),list(df.loc[df['sample'] == 1]["fold_change"]))[1], precision=1, unique=False, fractional=False, trim='k')
                axes[chr].plot([sample_number, sample_number, number, number], [height+i+i+0.1, height+i+i+0.35, height+i+i+0.35, height+i+i+0.1], lw=1.5, c='k')
                axes[chr].text((sample_number+number)*.5, height+i+i+0.35, statistic, ha='center', va='bottom',fontsize=13)
    #axes[x_pos,y_pos].set_xlim(-0.5,3.5)
    plt.tight_layout()


#axes[2].legend(loc=("lower right"))
plt.savefig("{}{}".format(outputdir,'/Isoseq_foldchange_117_152.pdf')  )

