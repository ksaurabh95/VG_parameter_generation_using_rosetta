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
# output format = [thetar, thetas, alpha (1/cm), n, Ksat (cm/day) ]

sand = 30 
silt = 30
clay = 40

depth = [10]
CI = [2.5,97.5] 

textural_data =  [  [sand,silt,clay] ]

model = SoilHydraulicModel(textural_data,depth,CI)
print(model.textures)  # ['sandy loam', 'loam', 'clay loam']
mean, stdev = model.get_parameters()
model.plot_mean_SWRC_HCC_curve() 
model.plot_SWRC_HCC_with_confidence_band() 




    



