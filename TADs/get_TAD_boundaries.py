
# Generate cooler files from hic
from hic2cool import hic2cool_convert
hic2cool_convert("inter_WT.hic", "OESO_103_WT_chr6_assmebly_inter.cool")
hic2cool_convert("inter_CT.hic", "OESO_103_CT_chr6_assmebly_inter.cool")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import statistics

import cooler
import cooltools.lib.plotting
from cooltools import insulation

#from packaging import version
#version.parse(cooltools.__version__):
#raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
#                         "please check your cooltools version and update to the latest")

import cooltools

resolution = 10000
clr = cooler.Cooler(f'OESO_103_{chr}_chr6_assmebly_inter.mcool::resolutions/{resolution}')
#alter the resolution values to make compatible for plotting - ONLY NEED TO DO THIS ONCE
with clr.open('r+') as f:
    f["bins"].create_dataset("KR_mult", data=1./f["bins/KR"][:], compression="gzip", compression_opts=6)


#windows = [15*resolution]
windows = [15*resolution,25*resolution,50*resolution,75*resolution,100*resolution]

chromosomes=["WT","CT"]
for win in windows:
    print(win)
    for chr in chromosomes:
    # define the cool file
        clr = cooler.Cooler(f'OESO_103_{chr}_chr6_assmebly_inter.mcool::resolutions/{resolution}')
        # load in the cool file
        insulation_table = insulation(clr, win,verbose=True,clr_weight_name="KR_mult",ignore_diags=2)
        insulation_table.to_csv('{}{}{}{}{}{}{}'.format("insulation_table_",chr,"_resolution_",resolution,"_windows_",win,".txt"),index=False, sep="\t")
        # read in the contig lengths file and turn it into a dictionary
        index=pd.read_csv("{}{}{}".format("/lustre/scratch126/casm/team154pc/ji2/OesophagealLongRead_WTSI-OESO_103/hifiasm_assembly/hifiasm_assembly/final_hifiasm_assembly/assembly/chr6-hifiasm-caus3D-",chr,".renamed.contiglengths.txt",sep=""), sep="\t", header=None)
        index=index.set_axis(['name', 'length'], axis=1)
        index['name']=index['name'].str.upper()
        index = index.set_index('name').T.to_dict('list') 
        # get the insulation table with the boundary that you want 
        insulation_table_filt= insulation_table[insulation_table["{}{}".format('is_boundary_',win)]].reset_index()
        insulation_table_filt["last_tad_size"]=0
        insulation_table_filt["next_tad_size"]=0
        for row in range(len(insulation_table_filt)):
            if row==0:
                insulation_table_filt.at[insulation_table_filt.index[row], 'last_tad_size']=insulation_table_filt["start"][row]
            elif insulation_table_filt["chrom"][row]==insulation_table_filt["chrom"][row-1]:
                insulation_table_filt.at[insulation_table_filt.index[row], 'last_tad_size']=insulation_table_filt["start"][row]-insulation_table_filt["start"][row-1]
            elif insulation_table_filt["chrom"][row]!=insulation_table_filt["chrom"][row-1]:
                insulation_table_filt.at[insulation_table_filt.index[row], 'last_tad_size']=insulation_table_filt["start"][row]
        for row in range(len(insulation_table_filt)):
            if row==len(insulation_table_filt)-1:
                insulation_table_filt.at[insulation_table_filt.index[row], 'next_tad_size']=index[insulation_table_filt["chrom"][row]]-insulation_table_filt["start"][row]
            elif insulation_table_filt["chrom"][row]==insulation_table_filt["chrom"][row+1]:
                insulation_table_filt.at[insulation_table_filt.index[row], 'next_tad_size']=insulation_table_filt["start"][row+1]-insulation_table_filt["start"][row]
            elif insulation_table_filt["chrom"][row]!=insulation_table_filt["chrom"][row+1]:
                insulation_table_filt.at[insulation_table_filt.index[row], 'next_tad_size']=index[insulation_table_filt["chrom"][row]]-insulation_table_filt["start"][row]
        insulation_table_filt.to_csv('{}{}{}{}{}{}{}'.format("insulation_table_",chr,"_resolution_",resolution,"_windows_",win,"_with_boundaries.txt"),index=False, sep="\t")
