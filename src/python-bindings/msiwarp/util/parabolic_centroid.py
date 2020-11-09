# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 09:41:33 2019

@author: jo5130er
"""

from scipy.signal import find_peaks
import numpy as np

def parabolic_centroid(mzs, intensities, peak_threshold):
    """Centroid
    
    returns:
    """    
    peak_indices, _ = find_peaks(intensities, height=peak_threshold)
    peak_left = peak_indices - 1
    peak_right = peak_indices + 1
    
    n = len(peak_indices)
    
    X = np.zeros((n, 3))
    Y = np.zeros((n, 3))
    
    X[:,0] = mzs[peak_left]
    X[:,1] = mzs[peak_indices]
    X[:,2] = mzs[peak_right]
    
    Y[:,0] = intensities[peak_left]
    Y[:,1] = intensities[peak_indices]
    Y[:,2] = intensities[peak_right]
    
    a = ((Y[:,2] - Y[:,1]) / (X[:,2] - X[:,1]) - 
         (Y[:,1] - Y[:,0]) / (X[:,1] - X[:,0])) / (X[:,2] - X[:,0])
    
    b = ((Y[:,2] - Y[:,1]) / (X[:,2] - X[:,1]) * (X[:,1] - X[:,0]) + 
         (Y[:,1] - Y[:,0]) / (X[:,1] - X[:,0]) * (X[:,2] - X[:,1])) / (X[:,2] - X[:,0])             

    mzs_parabolic = ((1/2) * (-b + 2 * a * X[:,1]) / a)
    intensities_parabolic = (a * (mzs_parabolic - X[:,1]) ** 2 +
                             b * (mzs_parabolic - X[:,1]) + Y[:,1])
    
    return (mzs_parabolic, intensities_parabolic)