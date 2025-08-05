# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 14:21:56 2025

@author: Saurabh
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rosetta import rosetta, SoilData



def classify_usda_texture(sand, silt, clay):
    if (sand + silt + clay) != 100:
        raise ValueError("Percentages must sum to 100.")

    if clay >= 40 and silt >= 40:
        return 'silty clay'
    elif clay >= 35 and sand >= 45:
        return 'sandy clay'
    elif clay >= 35:
        return 'clay'
    elif silt >= 80:
        return 'silt'
    elif silt >= 50 and clay < 27:
        return 'silt loam'
    elif sand >= 85:
        return 'sand'
    elif sand >= 70 and clay < 15:
        return 'loamy sand'
    elif sand >= 52 and clay < 20:
        return 'sandy loam'
    elif clay >= 27 and clay < 40 and sand > 20 and sand < 45:
        return 'clay loam'
    elif clay >= 20 and clay < 35 and silt > 27 and silt < 50:
        return 'loam'
    elif clay < 27 and silt >= 50 and silt < 80:
        return 'silty loam'
    else:
        return 'loam'  # Fallback
    

def Van_Genuchten_moisture(h, thetar, thetas, alpha, N):
    """
    Calculate soil moisture based on the Van Genuchten relationship (1981)
    h: pressure head (can be scalar or numpy array)
    thetar: residual moisture content
    thetas: saturated moisture content
    alpha, N: van Genuchten parameters
    """
    h = np.asarray(h)  # ensures compatibility with scalar and array input
    m = 1 - 1 / N

    # For h <= 0 (unsaturated zone)
    beta = np.power(np.abs(alpha * h), N)
    theta_unsat = thetar + (thetas - thetar) * np.power(1 + beta, -m)

    # Saturated condition where h > 0
    # theta = np.where(h >= 0, theta_unsat, thetas)
    theta = theta_unsat 

    return theta



def Van_Genuchten_K(h, Ksat, alpha, N):
    """
    Unsaturated hydraulic conductivity based on modified Van Genuchten–Mualem model (Nelson 1985)
    
    h: pressure head (can be scalar or numpy array)
    Ksat: saturated hydraulic conductivity
    thetar, thetas: residual and saturated moisture content
    alpha, N: van Genuchten parameters
    n_eta: relative permeability exponent , 0.5 as default
    """
    h = np.asarray(h)
    m = 1 - 1 / N

    # Compute Se
    beta = np.power(np.abs(alpha * h), N)
    Se = np.power(1 + beta, -m)

    # Compute K(h) using np.where for condition h < 0
    term1 = np.power((1 - np.power(1 - np.power(Se, 1/m), m)), 2)
    K_unsat = Ksat * np.power(Se, 0.5) * term1

    # Saturated where h >= 0
    K = np.where(h > 0, K_unsat, Ksat)

    return K



class SoilHydraulicModel:
        
    def __init__(self, textural_data, depth=None,CI =None ):
        self.depth = depth
        self.textural_data = textural_data
        self.CI = CI

        self.num_layers = len(self.textural_data)

        self.textures = []  # e.g., ['loam', 'silt loam', ...]
 
        for i, row in enumerate(self.textural_data):
            sand, silt, clay = row
            texture = classify_usda_texture(sand, silt, clay)
            self.textures.append(texture)
        
    def get_parameters(self):
        
        soildata = SoilData.from_array(self.textural_data)
        
        mean, stdev, codes = rosetta(3, soildata)
        mean = np.atleast_2d(mean)
        stdev = np.atleast_2d(stdev)

        # mean[:,2:6] = pow(10,mean[:,2:6] ) 
        # stdev[:,2:6] = pow(10,stdev[:,2:6] ) 
   
        self.params_mean_list = mean
        self.params_std_list = stdev
        
        param_names = ['thetar', 'thetas', 'log10(alpha) (1/cm)', 'log10(n)', 'log10(Ksat) (cm/day)']
        
        # parameters = [mean, wilting_point, field_capacity ]
        self.df_mean = pd.DataFrame(mean, columns=param_names)
        self.df_std = pd.DataFrame(stdev, columns=param_names)
        
        return self.df_mean, self.df_std
    
    def wilting_point_and_field_calacity(self):
        # wilting point and field capacity is estimated when pressure head is 330 cm and 15000 cm respectively
        h_wp = -15000 # pressure head value in cm at permanent wilting point 
        h_fc = -300  # pressure head value in cm at field capacity point
        
        mean = self.params_mean_list
        depth = self.depth
        mean = np.atleast_2d(mean)
        wilting_point, field_capacity = [], []
        
        
        for i in range(len(depth)):
            thetar, thetas, log10_alpha,log10_n, log10_Ksat = mean[i,:]
            alpha = pow(10,log10_alpha)
            n = pow(10,log10_n)
            # θ(h)
            theta_wp = Van_Genuchten_moisture(h_wp, thetar, thetas, alpha, n)
            theta_fc = Van_Genuchten_moisture(h_fc, thetar, thetas, alpha, n)
            
            wilting_point.append(theta_wp)
            field_capacity.append(theta_fc)
        
        # param_names = ['wilting_point', 'field capacity']
        # pramas = [wilting_point,field_capacity]    
        wilting_point = np.array(wilting_point)
        field_capacity = np.array(field_capacity)
        param_names = ['wilting_point', 'field capacity']
        pramas = np.array([wilting_point,field_capacity] ) 
        pramas = pramas.transpose()
        
        self.wp_fc = pd.DataFrame(pramas, columns=param_names)
        
        # self.wilting_point = pd.DataFrame(wilting_point, columns='wilting_point')
        # self.field_capacity = pd.DataFrame(field_capacity, columns='field_capacity')

        
         
        return self.wp_fc 
    
    
    
    def plot_mean_SWRC_HCC_curve(self):
        
        h = np.logspace(0, 5, 200)
        
        mean = self.params_mean_list
        depth = self.depth
        mean = np.atleast_2d(mean)
        texture = self.textures


        # HCC and SWRC curve will be plotted here 
        layers = len(depth) 
        base_width=5
        base_height=4
        width = base_width * 3
        height = base_height * layers
        fig, axs = plt.subplots( layers,3,figsize=(width, height), squeeze=False)
        # fig, axs = plt.subplots( layers,3,figsize=(8*layers,5))

        # Fix for subplot axis layout
        if layers == 1:
            axs = axs.reshape(1, -1)  # force to 2D shape (1, 3)


        for i in range(len(depth)):
            thetar, thetas, log10_alpha,log10_n, log10_Ksat = mean[i,:]
            alpha = pow(10,log10_alpha)
            n = pow(10,log10_n)
            Ksat = pow(10,log10_Ksat)

            # θ(h)
            theta = Van_Genuchten_moisture(h, thetar, thetas, alpha, n)
            # K(h)
            K =  Van_Genuchten_K(h, Ksat, alpha, n) 
        
            axs[i,0].plot(h,theta)
            axs[i,0].set_xscale('log')
            axs[i,0].set_xlabel('|h| (cm)')
            axs[i,0].set_ylabel(r'$\theta$')
            axs[i,0].set_title(f'θ(h) vs |h| Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,0].set_ylim([0,1])
            axs[i,0].set_xlim([1,10**5])
            axs[i,0].grid(True)
    
    
            axs[i,1].plot(h,K)
            axs[i,1].set_xscale('log')
            axs[i,1].set_xlabel('|h| (cm)')
            axs[i,1].set_ylabel(r'$K$ $(cm/day)$')
            axs[i,1].set_title(f'K(h) vs |h| Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,1].set_xlim([1,10**5])
            axs[i,1].set_ylim([0,np.round(Ksat)])

            axs[i,1].grid(True)
    
            axs[i,2].plot(theta,K)
            axs[i,2].set_xlabel(r'$\theta$')
            axs[i,2].set_ylabel(r'$K$ $(cm/day)$')
            axs[i,2].set_title(f'K(h) vs θ(h) Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,2].set_xlim([0,1])
            axs[i,2].set_ylim([0,np.round(Ksat)])

            axs[i,2].grid(True)

        plt.tight_layout()
        plt.show()
        
    def _monte_carlo_samples(self, h,mean ,stdev,N=1000):
        # monte carlo analysis will be performed for each soil layer 1000 times 
        
        samples = np.random.normal(loc=mean, scale=stdev, size=(N, 5))

        thetar_mean, thetas_mean, log10_alpha_mean,log10_n_mean, log10_Ksat_mean = mean
        alpha_mean = pow(10,log10_alpha_mean )
        n_mean = pow(10,log10_n_mean )
        Ksat_mean = pow(10, log10_Ksat_mean)
        # θ(h) mean
        theta_mean  = Van_Genuchten_moisture(h, thetar_mean , thetas_mean , alpha_mean , n_mean )
        # K(h) mean
        K_mean  =  Van_Genuchten_K(h, Ksat_mean , alpha_mean , n_mean )
               
        theta_samples, K_samples = [], []
        
        for j in range(N):
            thetar, thetas, log10_alpha,log10_n, log10_Ksat = samples[j,:]
            alpha = pow(10,log10_alpha )
            n = pow(10,log10_n )
            Ksat = pow(10, log10_Ksat)
            theta = Van_Genuchten_moisture(h, thetar, thetas, alpha, n)
            # K(h)
            K =  Van_Genuchten_K(h, Ksat, alpha, n)
            
            theta_samples.append(theta)
            K_samples.append(K)
            
        theta_samples = np.array(theta_samples) 
        K_samples = np.array(K_samples)
        
        return theta_mean,theta_samples,K_mean,K_samples
    
        
    
    def plot_SWRC_HCC_with_confidence_band(self):
        # CI refer to the confidence interval, [2.5,97.5] for 95 % CI
        CI = self.CI
        lower_ci = CI[0]
        upper_ci = CI[1]
        
        h = np.logspace(0, 5, 200)
        
        mean = self.params_mean_list
        stdev = self.params_std_list
        depth = self.depth
        texture = self.textures

        mean = np.atleast_2d(mean)
        stdev = np.atleast_2d(stdev)
        # HCC and SWRC curve will be plotted here 
        layers = len(depth) 
        base_width=5
        base_height=4
        width = base_width * 3
        height = base_height * layers
        fig, axs = plt.subplots( layers,3,figsize=(width, height), squeeze=False)
        
        if layers == 1:
            axs = axs.reshape(1, -1)  # force to 2D shape (1, 3)
        for i in range(len(depth)):

            theta_mean,theta_samples,K_mean,K_samples = self._monte_carlo_samples(h, mean[i,:],  stdev[i,:] ) 
            
            theta_lower = np.percentile(theta_samples, lower_ci, axis=0)
            theta_upper = np.percentile(theta_samples, upper_ci, axis=0)       
            
            K_lower = np.percentile(K_samples, lower_ci, axis=0)
            K_upper = np.percentile(K_samples, upper_ci, axis=0)
            
            
            axs[i,0].plot(h,theta_mean, label='mean')
            axs[i,0].fill_between(h, theta_lower, theta_upper, color='lightgreen', alpha=0.5, label='95% CI')
            axs[i,0].set_xscale('log')
            axs[i,0].set_xlabel('|h| (cm)')
            axs[i,0].set_ylabel(r'$\theta$')
            axs[i,0].set_title(r'$\theta$ vs |h| at depth')
            axs[i,0].set_title(f'θ(h) vs |h| Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,0].set_ylim([0,1])
            axs[i,0].set_xlim([1,10**5])
            axs[i,0].grid(True)
    
    
            axs[i,1].plot(h,K_mean, label='mean')
            axs[i,1].fill_between(h, K_lower, K_upper, color='lightgreen', alpha=0.5, label='95% CI')
            axs[i,1].set_xscale('log')
            axs[i,1].set_xlabel('|h| (cm)')
            axs[i,1].set_ylabel(r'$K$ $(cm/day)$')
            axs[i,1].set_title(f'K(h) vs |h| Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,1].set_xlim([1,10**5])

            axs[i,1].grid(True)
    
            axs[i,2].plot(theta_mean,K_mean, label='mean')
            axs[i,2].fill_between(theta_mean, K_lower, K_upper, color='lightgreen', alpha=0.5, label='95% CI')

            axs[i,2].set_xlabel(r'$\theta$')
            axs[i,2].set_ylabel(r'$K$ $(cm/day)$')
            axs[i,2].set_title(f'K(h) vs θ(h) Depth: {depth[i]} cm for {texture[i]} soil') 
            axs[i,2].set_xlim([0,1])

            axs[i,2].grid(True)

        plt.tight_layout()
        plt.show()


        
        
        
                
                    
                    
                    
                    
                
                
                
                
                
                
                
  



    
    
    
    
    
    
    
    
    
    
    
    
    
