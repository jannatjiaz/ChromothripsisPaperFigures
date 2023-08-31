
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from numpy import *

import cooler
import cooltools.lib.plotting
from cooltools import insulation

from packaging import version

from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe

# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

plt.rcParams['font.size'] = 12

chr="CT"
resolution = 10000
windows = [100*resolution]
windows = [15*resolution]

clr = cooler.Cooler(f'OESO_103_{chr}_chr6_assmebly_inter.mcool::resolutions/{resolution}')
insulation_table = insulation(clr, windows,verbose=True,clr_weight_name="KR_mult",ignore_diags=2)

index=pd.read_csv("{}{}{}".format("/lustre/scratch126/casm/team154pc/ji2/OesophagealLongRead_WTSI-OESO_103/hifiasm_assembly/hifiasm_assembly/final_hifiasm_assembly/assembly/chr6-hifiasm-caus3D-",chr,".renamed.contiglengths.txt",sep=""), sep="\t", header=None)

start = 50_000_000
end = start+ 7443939  #max is 25000000

for row in range(len(index)):
#for row in range(1,len(index)):
    chromosome_name=index.loc[row][0].upper()
    chromosome_name_lower=index.loc[row][0]
    start = 0
    end = start+ index.loc[row][1]  #max is 25000000
    region = (chromosome_name, start, end)
    norm = LogNorm(vmax=5, vmin=0.3) #change the colours of the hic matrix
    data = clr.matrix(balance=False).fetch(region)
    f, ax = plt.subplots(figsize=(20, 10))
    im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='YlOrRd')
    #im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
    ax.set_aspect(0.5)
    ax.set_ylim(0, 50*4.5*resolution) # change the height of the plots 
    format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0., aspect=10)
    plt.colorbar(im, cax=cax)
    insul_region = bioframe.select(insulation_table, region)
    SVs=pd.read_csv("{}{}{}".format("/lustre/scratch126/casm/team154pc/ji2/OesophagealLongRead_WTSI-OESO_103/hifiasm_assembly/hifiasm_assembly/final_hifiasm_assembly/index/",chr,"_hifiasm_assembly_chr6_all_contigs_index_for_TADs.txt"), sep="\t", header=None)
    SVs=SVs.set_axis(['assembly_chr', 'assembly_start', 'assembly_end','ref_chr', 'ref_start', 'ref_end'], axis=1)
    SVs=SVs.loc[SVs['assembly_chr'] == chromosome_name_lower]
    SVs=SVs.loc[SVs['assembly_start'] > start]
    SVs=SVs.loc[SVs['assembly_end'] < end]
    SVs=SVs.reset_index() 
    # line 1 points
    ins_ax = divider.append_axes("bottom", size="10%", pad=0., sharex=ax)
    for x in range(len(SVs)):
        if SVs.loc[x]["assembly_end"]>SVs.loc[x]["assembly_start"]:
            x1 = [SVs.loc[x]["assembly_start"],SVs.loc[x]["assembly_start"]]
            y1 = [0,1]
        else:
            x1 = [SVs.loc[x]["assembly_end"],SVs.loc[x]["assembly_end"]]
            y1 = [0,1]       
        # plotting the line 1 points 
        ins_ax.plot(x1, y1, label = "line 1", color="black")
    ins_ax.yaxis.set_visible(False)
    ins_ax.xaxis.set_visible(False)
    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
    #ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                insul_region['log2_insulation_score_'+str(windows[0])],
                label=f'Window {windows[0]} bp')
    boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
    weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
    strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]
    ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
                weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
    ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
                strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')
    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_xlim(region[1], region[2])
    plt.savefig('{}{}{}{}{}{}'.format(chr,'_boundaries_called_',windows[0],'_',chromosome_name,'.pdf'))
    #plt.show()

