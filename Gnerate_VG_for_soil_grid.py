# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 11:03:57 2025

@author: Saurabh
"""

import numpy as np 
from rosetta import rosetta, SoilData
import matplotlib.pyplot as plt
import scipy.io
from functions_required import SoilHydraulicModel

# input format = [%sand, %silt, %clay, buld density, th33, th1500] 
# input format = [%sand, %silt, %clay] 

# output format = [thetar, thetas, log10 alpha (1/cm), log10 (n), log10(Ksat) (cm/day) ]

sand = 30 
silt = 30
clay = 40

depth = [10,20] # depth in cm 
CI = [2.5,97.5] # confidence band or percentile 

textural_data =  [  [sand,silt,clay],[sand,silt,clay] ]

model = SoilHydraulicModel(textural_data,depth,CI)
print(model.textures)  # ['sandy loam', 'loam', 'clay loam']
mean_params, std_params = model.get_parameters()  # mean and standard deviation of vg parameters
# model.plot_mean_SWRC_HCC_curve() 
model.plot_SWRC_HCC_with_confidence_band() 

 







