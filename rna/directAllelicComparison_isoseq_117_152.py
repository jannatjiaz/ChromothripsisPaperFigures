
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math 
import numpy as np
from scipy.stats import pearsonr
import argparse 

parser = argparse.ArgumentParser(description='usage:  DirectAllelicComparison_isoseq_117_152.py   --outputdir outputdir')
parser.add_argument('--outputdir', help='path to output directory',required=True)
args = parser.parse_args()

outputdir=args.outputdir

chromosomes=["chr9","chr6"]
order=["WTSI-OESO_117","WTSI-OESO_152"]

color_dict = dict({'chromothriptic':'teal',
                  'wildtype':'darkmagenta'})


hap1_coverage_117_chr6=22771.7
hap1_coverage_152_chr6=20809.6
hap1_coverage_117_chr9=17557.7
hap1_coverage_152_chr9=17520.8
hap2_coverage_117_chr6=24465.7
hap2_coverage_152_chr6=27386.4
hap2_coverage_117_chr9=18590.4
hap2_coverage_152_chr9=21391.5

for chr in range(len(chromosomes)):
    chromosome=chromosomes[chr]
    chromosome
    #get input files paths
    input_counts=[]
    input_sites=[]
    for sample in order:
        s=sample.split("_")[1]
        input_sites.append("{}{}{}{}{}".format("example_inputs/differentialExpression_",chromosome,"_",s,"_filtered_final_recorrected.csv"))        
        input_counts.append("{}{}{}".format("example_inputs/",chromosomes[chr],"_counts_117_152.csv"))
    #add the differential gene names to the dictionary
    genes_list=[]
    for sample in range(len(order)):
        genes = pd.read_csv(input_sites[sample],sep="\t",header=None,skiprows=1)
        genes.columns = ["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr"]
        genes_list.append(list(genes["gene"]))
    genes_final = [value for value in genes_list[0] if value in genes_list[1]]
    genes_dict={x:{} for x in genes_final} #make the contig keys
    #Get the fold changes for the genes:    
    for sample in range(len(order)):
        sample
        genes = pd.read_csv(input_counts[sample],sep="\t",header=None,skiprows=1)
        #if differential == "differential":
        genes.columns = ["gene","117_1","117_2","152_1","152_2"]
        for item in range(len(genes)):
            if genes["gene"][item] in genes_dict:
                #normalise by differences in copy number 
                if chromosome=="chr6" and order[sample]=="WTSI-OESO_117":
                    ct = (genes["117_1"][item]/hap1_coverage_117_chr6)*hap1_coverage_152_chr6
                    wt = (genes["117_2"][item]/hap2_coverage_117_chr6)*hap2_coverage_117_chr6
                elif chromosome=="chr9" and order[sample]=="WTSI-OESO_117":
                    ct = (genes["117_1"][item]/hap1_coverage_117_chr9)*hap1_coverage_152_chr9
                    wt = (genes["117_2"][item]/hap2_coverage_117_chr9)*hap2_coverage_117_chr9
                else:
                    ct = genes["152_1"][item]
                    wt = genes["152_2"][item]
                if ct>1:
                    ct = math.log2(ct)
                else:
                    ct=0    
                if wt>1:
                    wt = math.log2(wt)  
                else:
                    wt=0   
                genes_dict[genes["gene"][item]][order[sample]]={"wt":wt,"ct":ct}
    x=[]
    y=[]
    col=[]
    for gene in genes_dict:
        if len(genes_dict[gene])==2:
            x.append(genes_dict[gene]['WTSI-OESO_117']['ct'])
            y.append(genes_dict[gene]['WTSI-OESO_152']['ct'])
            col.append("chromothriptic")
            x.append(genes_dict[gene]['WTSI-OESO_117']['wt'])
            y.append(genes_dict[gene]['WTSI-OESO_152']['wt'])
            col.append("wildtype")
    noise = np.random.normal(0,0.15,len(y))
    jittered_y = y + noise
    noise = np.random.normal(0,0.15,len(x))
    jittered_x = x + noise
    sns.set(font_scale = 2)
    sns.set_style(style='white') 
    g = sns.JointGrid(x=jittered_x,y=jittered_y,height=8,ratio=5,xlim=(-1, 12), ylim=(-1, 12))
    # add scatter plot layer
    g.plot_joint(sns.scatterplot,s=10,linewidth=0,hue=col,palette=color_dict)
    # add marginal density plot layer
    g.plot_marginals(sns.kdeplot,linewidth=5,color="black")
    g.fig.suptitle(chromosomes[chr], y=1.3)
    # add 0 line
    g.ax_joint.axvline(x=0, color='black', alpha=0.5)
    g.ax_joint.axhline(y=0, color='black', alpha=0.5)
    # JointGrid has a convenience function
    g.set_axis_labels('log2 normalised read count WTSI-OESO_117', 'log2 normalised read count WTSI-OESO_152')
    g.ax_joint.plot([0, 1], [0, 1], ':k', transform=g.ax_joint.transAxes)
    # labels appear outside of plot area, so auto-adjust
    #plt.tight_layout()
    #plt.yticks(np.arange(0, 20, 2.5))
    #plt.xticks(np.arange(0, 20, 2.5))
    plt.rcParams['pdf.fonttype'] = 42
    corr, _ = pearsonr(x, y)
    print(chromosomes[chr])
    print('Pearsons correlation: %.3f' % corr)
    plt.savefig("{}{}{}{}".format(outputdir,'/Isoseq_117_152_',chromosomes[chr],'_comparison.pdf')  )


