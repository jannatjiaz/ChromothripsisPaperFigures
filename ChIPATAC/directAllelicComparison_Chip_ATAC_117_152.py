
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math 
import numpy as np
from scipy.stats import pearsonr
import argparse 

parser = argparse.ArgumentParser(description='usage:  DirectAllelicComparison_117_152.py   --outputdir outputdir')
parser.add_argument('--outputdir', help='path to output directory',required=True)
args = parser.parse_args()

outputdir=args.outputdir

whole_genome_coverage_117=(8.3+7.0)/2
whole_genome_coverage_152=(7.3+6.3)/2

chromosomes=["chr9","chr6"]
order=["WTSI-OESO_117","WTSI-OESO_152"]

color_dict = dict({'chromothriptic':'teal',
                  'wildtype':'darkmagenta'})

for chr in range(len(chromosomes)):
    chromosome=chromosomes[chr]
    chromosome
    #get input files paths
    input_counts=[]
    input_sites=[]
    for sample in order:
        input_counts.append("{}{}{}{}{}".format("example_inputs/",sample,"_",chromosomes[chr],"_samtools_counts.txt"))
        input_sites.append("{}{}{}{}{}".format("example_inputs/",sample,"_",chromosomes[chr],"_filtered_final_recorrected.csv"))
    #add the differential gene names to the dictionary
    for sample in range(len(order)):
        peaks = pd.read_csv(input_sites[sample],sep="\t",header=None,skiprows=1)
        peaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr"]
        peak_ids=list(peaks["chr"].astype(str)+"_"+peaks["start"].astype(str)+"_"+peaks["end"].astype(str))
    peaks_dict={x:{} for x in peak_ids} #make the contig keys
    #Get the fold changes for the genes:    
    for sample in range(len(order)):
        foldchange = pd.read_csv(input_counts[sample],sep="\t",header=None,skiprows=1)
        #if differential == "differential":
        foldchange.columns = ["chr","start","end","hap1_CTCF_1","hap1_H3K27me3_1","hap1_H3K27ac_1","hap1_H3K4me3_1","hap1_CTCF_2","hap1_H3K27me3_2","hap1_H3K27ac_2","hap1_H3K4me3_2","hap1_F1_2","hap1_F2_2","hap2_CTCF_1","hap2_H3K27me3_1","hap2_H3K27ac_1","hap2_H3K4me3_1","hap2_CTCF_2","hap2_H3K27me3_2","hap2_H3K27ac_2","hap2_H3K4me3_2","hap2_F1_2","hap2_F2_2"]
        foldchange_ids=list(foldchange["chr"].astype(str)+"_"+foldchange["start"].astype(str)+"_"+foldchange["end"].astype(str))
        for item in range(len(foldchange_ids)):
            for x in peaks_dict.keys():
                if foldchange_ids[item].split("_")[0]==x.split("_")[0] and int(foldchange_ids[item].split("_")[-2])-1000 < int(x.split("_")[-2]) < int(foldchange_ids[item].split("_")[-2])+1000 and int(foldchange_ids[item].split("_")[-1])-1000 < int(x.split("_")[-1]) < int(foldchange_ids[item].split("_")[-1])+1000:
                    counts1=(foldchange['hap1_CTCF_1'][item]+foldchange['hap1_H3K27me3_1'][item]+foldchange['hap1_H3K27ac_1'][item]+foldchange['hap1_H3K4me3_1'][item]+foldchange['hap1_CTCF_2'][item]+foldchange['hap1_H3K27me3_2'][item]+foldchange['hap1_H3K27ac_2'][item]+foldchange['hap1_H3K4me3_2'][item]+foldchange['hap1_F1_2'][item]+foldchange['hap1_F2_2'][item])/10
                    counts2=(foldchange['hap2_CTCF_1'][item]+foldchange['hap2_H3K27me3_1'][item]+foldchange['hap2_H3K27ac_1'][item]+foldchange['hap2_H3K4me3_1'][item]+foldchange['hap2_CTCF_2'][item]+foldchange['hap2_H3K27me3_2'][item]+foldchange['hap2_H3K27ac_2'][item]+foldchange['hap2_H3K4me3_2'][item]+foldchange['hap2_F1_2'][item]+foldchange['hap2_F2_2'][item])/10
                    #correct for whole chromosome gains
                    if chromosome=="chr9" and order[sample]=="WTSI-OESO_117":
                        hap1_counts = (counts1)
                        hap2_counts = (counts2/2)
                    elif chromosome=="chr6" and order[sample]=="WTSI-OESO_117":
                        hap1_counts = (counts1)
                        hap2_counts = (counts2)
                    else:
                        hap1_counts = (counts1/whole_genome_coverage_152)*whole_genome_coverage_117
                        hap2_counts = (counts2/whole_genome_coverage_152)*whole_genome_coverage_117
                    if hap1_counts!=0:
                        hap1_counts = math.log2(hap1_counts)
                    if hap2_counts!=0:
                        hap2_counts = math.log2(hap2_counts)   
                    peaks_dict[x][order[sample]]={"hap1":hap1_counts,"hap2":hap2_counts}
    x=[]
    y=[]
    col=[]
    for peak in peaks_dict:
        if len(peaks_dict[peak])==2:
            x.append(peaks_dict[peak]['WTSI-OESO_117']['hap1'])
            y.append(peaks_dict[peak]['WTSI-OESO_152']['hap1'])
            col.append("chromothriptic")
            x.append(peaks_dict[peak]['WTSI-OESO_117']['hap2'])
            y.append(peaks_dict[peak]['WTSI-OESO_152']['hap2'])
            col.append("wildtype")
    #jittered_y = y + 0.5 * np.random.rand(len(y)) -0.05
    #jittered_x = x + 0.5 * np.random.rand(len(x)) -0.05    
    noise = np.random.normal(0,0.15,len(y))
    jittered_y = y + noise
    noise = np.random.normal(0,0.15,len(x))
    jittered_x = x + noise
    sns.set(font_scale = 2)
    sns.set_style(style='white') 
    g = sns.JointGrid(x=jittered_x,y=jittered_y,height=8,ratio=5,xlim=(-5, 12.5), ylim=(-5, 12.5))
    # add scatter plot layer
    g.plot_joint(sns.scatterplot,s=10,linewidth=0,hue=col,palette=color_dict)
    # add marginal density plot layer
    g.plot_marginals(sns.kdeplot,linewidth=5,color="black")
    g.fig.suptitle(chromosomes[chr], y=1.3)
    # add 0 line
    g.ax_joint.axvline(x=0, color='black', alpha=0.5)
    g.ax_joint.axhline(y=0, color='black', alpha=0.5)
    # JointGrid has a convenience function
    g.ax_joint.plot([0,1], [0,1], ':k', transform=g.ax_joint.transAxes)
    g.set_axis_labels('log2 raw read counts of WTSI-OESO_117', 'log2 raw read counts of WTSI-OESO_152')
    # labels appear outside of plot area, so auto-adjust
    plt.tight_layout()
    plt.rcParams['pdf.fonttype'] = 42
    corr, _ = pearsonr(x, y)
    print(chromosomes[chr])
    print('Pearsons correlation: %.3f' % corr)
    plt.savefig("{}{}{}{}".format(outputdir,'/Diffbind_117_152_',chromosomes[chr],'_comparison.pdf')  )

