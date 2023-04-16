
import pandas as pd
import argparse 
import math



parser = argparse.ArgumentParser(description='usage:  distanceofgenestopeaks.py  --differential differential --atacchip_diff atacchip_diff --atacchip_nondiff atacchip_nondiff --genes_input_dir genes_input_dir --genes_input_dir genes_input_dir --genes_file_name genes_file_name --ensemble_input ensemble_input')
parser.add_argument('--differential', help='differential or nonDifferential',required=True)
parser.add_argument('--atacchip_diff', help='path to diffferential peaks called from atac and chip-seq reads',required=True)
parser.add_argument('--atacchip_nondiff', help='path to non-diffferential peaks called from atac and chip-seq reads',required=True)
parser.add_argument('--genes_input_dir', help='path to directory containing differential and non-differentail genes',required=True)
parser.add_argument('--genes_file_name ', help='file name for differential or non differential peaks, must match differential perameter',required=True)
parser.add_argument('--ensemble_input ', help='path to ensemble gene ID file',required=True)
args = parser.parse_args()

differential=args.differential
atacchip_diff=args.atacchip_diff
atacchip_nondiff=args.atacchip_nondiff
genes_input_dir=args.genes_input_dir
genes_file_name=args.genes_file_name
ensemble_input=args.ensemble_input

genes_input = "{}{}{}".format(genes_input_dir,"/",genes_file_name)



#load sv calls - since this is huge keep it as a pandas df
diff_peak = pd.read_csv(atacchip_diff,sep="\t",header=None)
diff_peak.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","distance"]
diff_peak['center']=list(round(diff_peak["start"]+( (diff_peak["end"]- diff_peak["start"])/2 )))
#subset the sv calls for the lower coordinate
diff_peak_center = diff_peak[["chr","center"]]


#load sv calls - since this is huge keep it as a pandas df
nondiff_peak = pd.read_csv(atacchip_nondiff,sep="\t",header=None)
nondiff_peak.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","distance1","distance2"]
nondiff_peak['center']=list(round(nondiff_peak["start"]+( (nondiff_peak["end"]- nondiff_peak["start"])/2 )))
#subset the sv calls for the lower coordinate
nondiff_peak_center = nondiff_peak[["chr","center"]]

# load in the ensemble gene names 
ensemble_gene_info = pd.read_csv(ensemble_input,sep="\t",header=None, skiprows=1)
ensemble_gene_info_dict={x:{} for x in ensemble_gene_info[1]} #make the contig keys
for item in range(len(ensemble_gene_info)):
    chromosome = "chr"+str(ensemble_gene_info[2][item])
    if ensemble_gene_info[1][item]!=ensemble_gene_info[1][item]:
        ensemble_gene_info_dict[ensemble_gene_info[0][item]]={"chr":chromosome, "start":ensemble_gene_info[3][item], "end":ensemble_gene_info[4][item]}
    else:
        ensemble_gene_info_dict[ensemble_gene_info[1][item]]={"chr":chromosome, "start":ensemble_gene_info[3][item], "end":ensemble_gene_info[4][item]}


#load peak calls - since this is huge keep it as a pandas df
peaks = pd.read_csv(genes_input,sep="\t",header=None,skiprows=0)
#if differential == "differential":
peaks.columns =["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr"]

#turn peaks calls into a dictionary 
peaks_dict={x:{} for x in peaks["gene"]} #make the contig keys
for item in range(len(peaks_dict)):
    if  ensemble_gene_info_dict[peaks["gene"][item]]["chr"] == 'chr6':
        peaks_dict[peaks["gene"][item]]={"gene":peaks["gene"][item], "baseMean":peaks["baseMean"][item], "log2FoldChange":peaks["log2FoldChange"][item], "lfcSE":peaks["lfcSE"][item],"stat":peaks["stat"][item],"pvalue":peaks["pvalue"][item],"padj":peaks["padj"][item],"new_fdr":peaks["new_fdr"][item],"chr":ensemble_gene_info_dict[peaks["gene"][item]]["chr"],"start":ensemble_gene_info_dict[peaks["gene"][item]]["start"],"end":ensemble_gene_info_dict[peaks["gene"][item]]["end"]}

peaks_dict= {k: v for k, v in peaks_dict.items() if v}


#convert the start stop postions to reference positions
for peak in peaks_dict:
    peak
    peak_middle = (peaks_dict[peak]["start"]+peaks_dict[peak]["end"])/2
    chr = peaks_dict[peak]["chr"]
    tmp_breakpoints= diff_peak_center.loc[diff_peak_center['chr'] == chr]
    match = tmp_breakpoints.iloc[ (tmp_breakpoints['center']-peak_middle).abs().argsort()[:1] ]
    peaks_dict[peak]['distance_to_diffpeak_middle']=abs(peak_middle-int(match['center']))
    tmp_breakpoints= nondiff_peak_center.loc[nondiff_peak_center['chr'] == chr]
    match = tmp_breakpoints.iloc[ (tmp_breakpoints['center']-peak_middle).abs().argsort()[:1] ]
    peaks_dict[peak]['distance_to_nondiffpeak_middle']=abs(peak_middle-int(match['center']))

order=["chr","start","end","gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","new_fdr","distance_to_diffpeak_middle","distance_to_nondiffpeak_middle"]
output_peaks_df = pd.DataFrame(peaks_dict).T.reindex(columns=order)
output_peaks_df.to_csv("{}{}{}".format('nearestSV_',differential,'Peaks.csv'), index=False, sep='\t', na_rep="none", header=order)

