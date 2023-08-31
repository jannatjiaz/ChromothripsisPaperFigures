
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
