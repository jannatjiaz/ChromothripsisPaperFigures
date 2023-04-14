### nearest SV on the same allele

import pandas as pd
import argparse 
import math

parser = argparse.ArgumentParser(description='usage:  nearestSV.py  --differential differential ')
parser.add_argument('--differential', help='differential or nonDifferential',required=True)
parser.add_argument('--wt_sv_input', help='path to wild-type SVs file',required=True)
parser.add_argument('--Ct_sv_input', help='path to chromothriptic SVs file',required=True)
args = parser.parse_args()

differential=args.differential
wt_sv_input=args.wt_sv_input
wt_sv_input=args.ct_sv_input
peaks_input = "{}{}".format(differential,"Peaks.csv")



#load sv calls - since this is huge keep it as a pandas df
wt_sv = pd.read_csv(wt_sv_input,sep="\t",header=None)
wt_sv.columns = ["chr1","pos1","pos2","chr2","pos3","pos4"]
#wt_sv.columns = ["chr1","pos1","pos2","chr2","pos3","pos4","id","qual","strand1","strand2","ref","type","filter","info","format","info2"]
#subset the sv calls for the lower coordinate 
wt_sv_lower = wt_sv[["chr1","pos2"]]
#subset the sv calls for the upper coordinate 
wt_sv_upper = wt_sv[["chr2","pos4"]]
wt_sv_upper.columns=["chr1","pos2"]
# Stack the DataFrames on top of each other
wt_sv_breakpoints = pd.concat([wt_sv_lower, wt_sv_upper], axis=0)
# Reset the index values to the second dataframe appends properly
wt_sv_breakpoints = wt_sv_breakpoints.reset_index(drop=True)


ct_sv = pd.read_csv(ct_sv_input,sep="\t",header=None)
ct_sv.columns = ["chr1","pos1","pos2","chr2","pos3","pos4"]
#ct_sv.columns = ["chr1","pos1","pos2","chr2","pos3","pos4","id","qual","strand1","strand2","ref","type","filter","info","format","info2"]
#subset the sv calls for the lower coordinate 
ct_sv_lower = ct_sv[["chr1","pos2"]]
#subset the sv calls for the upper coordinate 
ct_sv_upper = ct_sv[["chr2","pos4"]]
ct_sv_upper.columns=["chr1","pos2"]
# Stack the DataFrames on top of each other
ct_sv_breakpoints = pd.concat([ct_sv_lower, ct_sv_upper], axis=0)
# Reset the index values to the second dataframe appends properly
ct_sv_breakpoints = ct_sv_breakpoints.reset_index(drop=True)


#load peak calls - since this is huge keep it as a pandas df
peaks = pd.read_csv(peaks_input,sep="\t",header=None,skiprows=0)
peaks.columns = ["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr"]
#turn peaks calls into a dictionary 
peaks_ids = peaks["chr"].astype(str)+"_"+peaks["start"].astype(str)
peaks_dict={x:{} for x in peaks_ids} #make the contig keys
for item in range(len(peaks_dict)):
    peaks_dict[peaks_ids[item]]={"chr":peaks["chr"][item], "start":peaks["start"][item], "end":peaks["end"][item], "width":peaks["width"][item],"strand":peaks["strand"][item],"Conc":peaks["Conc"][item],"Conc_chromothr":peaks["Conc_chromothr"][item],"Conc_wildtype":peaks["Conc_wildtype"][item],"Fold":peaks["Fold"][item],"p.value":peaks["p.value"][item],"FDR":peaks["FDR"][item],"new_fdr":peaks["new_fdr"][item]}


#convert the start stop postions to reference positions
for peak in peaks_dict:
    peak
    peak_middle = (peaks_dict[peak]["start"]+peaks_dict[peak]["end"])/2
    chr = peaks_dict[peak]["chr"]
    if differential=="differential":
        if peaks_dict[peak]["Conc_wildtype"]>peaks_dict[peak]["Conc_chromothr"]:
            tmp_breakpoints= wt_sv_breakpoints.loc[wt_sv_breakpoints['chr1'] == chr]
            match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
            peaks_dict[peak]['SV_distance_to_peak_middle']=abs(peak_middle-int(match['pos2']))
        else:
            tmp_breakpoints= ct_sv_breakpoints.loc[ct_sv_breakpoints['chr1'] == chr]
            match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
            peaks_dict[peak]['SV_distance_to_peak_middle']=abs(peak_middle-int(match['pos2']))
    else:
        tmp_breakpoints= wt_sv_breakpoints.loc[wt_sv_breakpoints['chr1'] == chr]
        match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
        peaks_dict[peak]['SV_distance_to_peak_middle_wt']=abs(peak_middle-int(match['pos2']))
        tmp_breakpoints= ct_sv_breakpoints.loc[ct_sv_breakpoints['chr1'] == chr]
        match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
        peaks_dict[peak]['SV_distance_to_peak_middle_ct']=abs(peak_middle-int(match['pos2']))



if differential == "differential":
    order=["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle"]
    output_peaks_df = pd.DataFrame(peaks_dict).T.reindex(columns=order)
elif differential == "nonDifferential":
    order=["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct"]
    output_peaks_df = pd.DataFrame(peaks_dict).T.reindex(columns=order)

output_peaks_df.to_csv("{}{}{}".format('nearestSV_',differential,'Peaks.csv'), index=False, sep='\t', na_rep="none", header=order)



### nearest SV on the opposite allele


#convert the start stop postions to reference positions
for peak in peaks_dict:
    peak
    peak_middle = (peaks_dict[peak]["start"]+peaks_dict[peak]["end"])/2
    chr = peaks_dict[peak]["chr"]
    if differential=="differential":
        if peaks_dict[peak]["Conc_wildtype"]<peaks_dict[peak]["Conc_chromothr"]:
            tmp_breakpoints= wt_sv_breakpoints.loc[wt_sv_breakpoints['chr1'] == chr]
            match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
            peaks_dict[peak]['SV_distance_to_peak_middle']=abs(peak_middle-int(match['pos2']))
        else:
            tmp_breakpoints= ct_sv_breakpoints.loc[ct_sv_breakpoints['chr1'] == chr]
            match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
            peaks_dict[peak]['SV_distance_to_peak_middle']=abs(peak_middle-int(match['pos2']))
    else:
        tmp_breakpoints= wt_sv_breakpoints.loc[wt_sv_breakpoints['chr1'] == chr]
        match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
        peaks_dict[peak]['SV_distance_to_peak_middle_ct']=abs(peak_middle-int(match['pos2']))
        tmp_breakpoints= ct_sv_breakpoints.loc[ct_sv_breakpoints['chr1'] == chr]
        match = tmp_breakpoints.iloc[ (tmp_breakpoints['pos2']-peak_middle).abs().argsort()[:1] ]
        peaks_dict[peak]['SV_distance_to_peak_middle_wt']=abs(peak_middle-int(match['pos2']))

if differential == "differential":
    order=["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle"]
    output_peaks_df = pd.DataFrame(peaks_dict).T.reindex(columns=order)
elif differential == "nonDifferential":
    order=["chr","start","end","width","strand","Conc","Conc_chromothr","Conc_wildtype","Fold","p.value","FDR","new_fdr","SV_distance_to_peak_middle_wt","SV_distance_to_peak_middle_ct"]
    output_peaks_df = pd.DataFrame(peaks_dict).T.reindex(columns=order)

output_peaks_df.to_csv("{}{}{}".format('nearestSV_',differential,'_on_other_allele_Peaks.csv'), index=False, sep='\t', na_rep="none", header=order)

