# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:03:55 2019

@author: jo5130er
"""

import numpy as np
import struct
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import savgol_filter

from .parabolic_centroid import parabolic_centroid


def update_meta(old_meta, spectra):
    """Generate new meta data for spectra with different number of points and
    offsets.
    
    Returns:
    list:List of meta objects
    """    
    offset = old_meta[0][0] # size of header and meta remains the same
    
    new_meta = old_meta.copy()    
    for i in range(len(old_meta)):
        n_points = len(spectra[i][0])
        
        m_i = old_meta[i] 
        new_meta[i] = (offset, n_points, m_i[2], m_i[3], m_i[4])
        
        offset = offset + n_points * 12
        
    return new_meta
    

def read_spectrum(filestream, n):
    """Reads n mz and intensity pairs from existing .sbd filestream. Assumes double
    precision for mass and single for intensity.
    
    Returns:
    Tuple:tuple of mz and intensity vectors    
    """    
    mzs = struct.unpack("<%dd" % n, filestream.read(8 * n))
    intensities = struct.unpack("<%df" % n, filestream.read(4 * n))    
    
    return (mzs, intensities)


def read_spectrum_fs(filepath, meta_item):
    """Reads n mz and intensity pairs from .sbd file stream. Assumes double
    precision for mass and single for intensity.
    
    Returns:
    Tuple:tuple of mz and intensity vectors 
    """
    with open(filepath, 'rb') as filestream:
        filestream.seek(meta_item[0])
        (mz, intensity) = read_spectrum(filestream, meta_item[1])
        
        return (np.array(mz), np.array(intensity))


from array import array
def write_spectrum(filestream, spectrum):
    """Writes n mz and intensity pairs from .sbd file stream. Assumes double
    precision for mass and single for intensity.    
    """    
    mzs = array('d', spectrum[0])
    intensities = array('f', spectrum[1])
    
    mzs.tofile(filestream)
    intensities.tofile(filestream)
    

def read_sbd_meta(filepath):
    with open(filepath, 'rb') as in_file:
        header = struct.unpack("<BQB", in_file.read(10))        
        meta_size = header[1] * 20 # sizeof(QLfHH)        
        meta = [meta_item for meta_item in 
                struct.iter_unpack("<QLfHH", in_file.read(meta_size))]
        return meta
    
    
def read_sbd(filepath):
    """Reads an .sbd file containing spectra in either profile or centroid mode    
    Returns:
    list:List of spectra    
    """    
    with open(filepath, 'rb') as in_file:
        header = struct.unpack("<BQB", in_file.read(10))
        
        meta_size = header[1] * 20 # sizeof(QLfHH)        
        meta = [meta_item for meta_item in
                struct.iter_unpack("<QLfHH", in_file.read(meta_size))]
        
        num_points = [meta_item[1] for meta_item in meta]        
        spectra = [read_spectrum(in_file, n) for n in num_points]
        
    return (header, meta, spectra)


def centroid_tof(spectrum, peak_threshold, window_size = 11, order = 2):
    """Centroids TOF spectrum by two rounds of Savitzy-Golay smoothing.    
    Returns:
    list: mz and intensity values of peak centroids
    """
    intensity_golay = savgol_filter(np.array(spectrum[1]), window_size, order)
    intensity_golay2 = savgol_filter(intensity_golay, 11, order)
    
    # Centroid smoothed signal
    (mz_c, intensity_c) = parabolic_centroid(np.array(spectrum[0]),
                                             intensity_golay2, peak_threshold)
    
    return (mz_c, intensity_c)    
      
    
def centroid_orbitrap(spectrum, peak_threshold):        
    """Centroids oribtrap spectrum and tries to remove artifacts from 
    spectral leakage.
    
    Returns:
    list: mz and intensity values of peak centroids
    """    
    mz_pks, pks = parabolic_centroid(np.array(spectrum[0]),
                                     np.array(spectrum[1]), peak_threshold)
    
    return (mzs_cluster, intensities_cluster)


def centroid_sbd(filepath_in, filepath_out, centroid_function, peak_threshold):
    """Reads an .sbd file containing spectra in either profile or centroid mode
    
    Returns:
    list:List of spectra    
    """    
    with open(filepath_in, 'rb') as in_file:
        header = struct.unpack("<BQB", in_file.read(10))
        
        meta_size = header[1] * 20 # sizeof(QLfHH)        
        meta = [meta_item for meta_item in
                struct.iter_unpack("<QLfHH", in_file.read(meta_size))]
        
        n_points = [meta_item[1] for meta_item in meta]
        
        centroided_spectra = []        
        for i in range(len(meta)):
            spectrum = read_spectrum(in_file, n_points[i])
            centroided_spectra.append(centroid_function(spectrum, peak_threshold))
            
        new_meta = update_meta(meta, centroided_spectra)
            
        # Write header, updated meta and centroided spectra to outfile
        with open(filepath_out, 'wb') as out_file:
            # Write data to stream...
            out_file.write(struct.pack('<BQB', header[0], header[1], header[2]))
            
            for meta_item in new_meta:
                out_file.write(struct.pack('<QLfHH',
                                           meta_item[0], meta_item[1], 
                                           meta_item[2], meta_item[3],
                                           meta_item[4]))
                
            for spectrum in centroided_spectra:
                write_spectrum(out_file, spectrum) 
           
            
class spectrum:
    def __init__(self, mzs, intensities, x, y, z):
        self.mzs = mzs
        self.intensities = intensities
        self.x = x
        self.y = y
        self.z = z

        
def import_imzml_dataset(filepath):
    """Reads an .imzml and stores    
    Returns:
    list:List of spectra    
    """
    p = ImzMLParser(filepath)
    
    spectra = []
    
    for idx, (x,y,z) in enumerate(p.coordinates):
        mzs, intensities = p.getspectrum(idx)
        spectra.append(spectrum(mzs, intensities, x, y, z))
        
    return spectra


def imzml_to_sbd(filepath_imzml, filepath_sbd):
    """Converts a pair of .imzml and .ibd files to .sbd   
      Returns:
      list:True on success    
      """    
    with open(filepath_sbd, 'wb') as out_file:
        p = ImzMLParser(filepath_imzml)
        n_spectra = len(p.coordinates)
        
        # First pass
        meta = []
        offset = 20 * n_spectra + 10       
        for idx, (x,y,z) in enumerate(p.coordinates):
            (mzs, intensities) = p.getspectrum(idx)
            n_points = len(mzs)
            
            meta.append((offset, n_points, np.sum(intensities), x, y))
            offset = offset + n_points * 12
        
        # Write data to stream...
        header = (0, n_spectra, 8)    
        out_file.write(struct.pack('<BQB', header[0], header[1], header[2]))    
        
        for meta_item in meta:
            out_file.write(struct.pack('<QLfHH',
                                       meta_item[0], meta_item[1], 
                                       meta_item[2], meta_item[3],
                                       meta_item[4]))
         
        # Second pass    
        for i in range(n_spectra):
            mzs, intensities = p.getspectrum(i)
            write_spectrum(out_file, (mzs, intensities)) 
    
    return True


def write_sbd(filepath, spectra, meta):
    with open(filepath, 'wb') as out_file:
        n_spectra = len(spectra)
        
        meta_new = []
        
        # First pass, update meta
        offset = 20 * n_spectra + 10 # Magic numbers       
        for idx, m_i in enumerate(meta):
            n_points = len(spectra[idx][0])
            (_, _, tic, x, y) = m_i
            
            meta_new.append((offset, n_points, tic, x, y))
            offset = offset + n_points * 12
            
        # Write data to stream...
        header = (0, n_spectra, 8)    
        out_file.write(struct.pack('<BQB', header[0], header[1], header[2]))
        
        for meta_item in meta_new:
            out_file.write(struct.pack('<QLfHH',
                                       meta_item[0], meta_item[1], 
                                       meta_item[2], meta_item[3],
                                       meta_item[4]))
            
        # Second pass    
        for (mzs, intensities) in spectra:
            write_spectrum(out_file, (mzs, intensities))
        
    return True