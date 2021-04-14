import numpy as np
import msiwarp as mx
import copy
import itertools
import matplotlib.pyplot as plt
from bisect import bisect_left
from scipy.signal import find_peaks
from .parabolic_centroid import parabolic_centroid
from .read_sbd import read_spectrum_fs

# ------  utility functions for the python interface ------
def get_mx_spectrum(file, meta, i, sigma_1, instrument_type):
    """Read a spectrum from .sbd file and set peak width according to sigma
    and instrument type.
    """
    s_i = read_spectrum_fs(file, meta[i])
    return to_mx_peaks(*s_i, sigma_1, id=i, instrument_type=instrument_type)


def to_mx_peaks(mzs, hs, sigma_1, id = 0, instrument_type = 'orbitrap'):
    """Generate a list of MSIWarp spectra, [s0, ..., sn], from a list of
     m/z values, mzs: [mz0, ..., mzn], and intensities, hs: [h0, ..., hn]. sigma_1 is
     the expected peak width at m/z = 1.0.
    """
    if(instrument_type == 'quadrupole'):
        return np.array([mx.peak(id, mz_i, intensity_i, sigma_1)
                     for (mz_i, intensity_i) in zip(mzs, hs)], dtype=mx.peak)
    elif(instrument_type == 'tof'):
        return  np.array([mx.peak(id, mz_i, intensity_i, sigma_1 * mz_i)
                     for (mz_i, intensity_i) in zip(mzs, hs)], dtype=mx.peak)
    elif(instrument_type == 'orbitrap'):
        return np.array([mx.peak(id, mz_i, intensity_i, sigma_1 * (mz_i ** (3/2)))
                     for (mz_i, intensity_i) in zip(mzs, hs)], dtype=mx.peak)
    elif(instrument_type =='ft-icr'):
        return  np.array([mx.peak(id, mz_i, intensity_i, sigma_1 * (mz_i ** 2))
                     for (mz_i, intensity_i) in zip(mzs, hs)], dtype=mx.peak)

def to_mz(mx_peaks):
    return np.array([p.mz for p in mx_peaks])

def to_height(mx_peaks):
    return np.array([p.height for p in mx_peaks])

def to_ion_image(s, meta):
    x = np.array([m_i[3] for m_i in meta])
    y = np.array([m_i[4] for m_i in meta])
    
    n = np.max(x)
    m = np.max(y)
    
    im = np.zeros((m, n))
    im[y - 1, x - 1] = s
    
    return im

def range_to_image(peaks, meta):
    s = np.zeros(len(meta))
    tic = np.array([m[2] for m in meta])    
    
    for p in peaks:
        index = p.id
        s[index] += p.height    
    
    return to_ion_image(s / tic, meta)


def dispersion_ppm(x):
    x_m = np.mean(x)
    x_sd = np.std(x)
    return x_sd / x_m * 1e6

def quantile_triplets(ts, q):
    mzs = np.array([t.mz for t in ts])
    heights = np.array([t.height for t in ts]) 
    sel = heights > np.quantile(heights, q)    
    return mzs[sel]

def dispersion_triplets(ts, q):
    mzs = quantile_triplets(ts, q)
    return dispersion_ppm(mzs)


def dispersion_sigma(x, sigma_0, instrument_type):
    x_m = np.mean(x)
    x_sd = np.std(x)
    
    if instrument_type == 'tof':
        return x_sd / (sigma_0 * x_m)
    elif instrument_type == 'orbitrap':
        return x_sd / (sigma_0 * x_m ** (3/2))

    
def quantile_peaks(peaks, q):
    mzs = to_mz(peaks)    
    heights = to_height(peaks)
    sel = heights > np.quantile(heights, q)    
    return mzs[sel]


def peak_density_mz(spectra, xi, bandwidth=15, threshold=0.1, stride=100):
    """Generate a density estimate of peak density on a set of m/z sampling points, xi,
    from a list of spectra. Returns the density curve and its peak x and y positions. 
    """
    yi = np.zeros(len(xi))

    for i in range(0, len(spectra), stride):
        s_i = spectra[i]
        p_i = [mx.peak(0, p.mz, 1, bandwidth) for p in s_i]
        yi += np.array(mx.splat_peaks(p_i, xi, 4.0))
        
    pks, _ = find_peaks(yi)   
    
    xp = xi[pks]
    yp = yi[pks]
    
    sel = yp > threshold * np.max(yp)
    
    return (yi, xp[sel], yp[sel])   


def generate_mean_spectrum(spectra, n_points, sigma_1, mz_begin,
                           mz_end, tic, instrument_type, stride = 100):
    """
    """
    xi_mean = np.zeros(n_points)
    xi_mean[0] = mz_begin
    
    d = 0.25 * sigma_1
    
    if instrument_type == 'tof':
        func = lambda x: x
    elif instrument_type == 'orbitrap':
        func = lambda x: x ** (3/2)
    elif instrument_type == 'ft-icr':
        func = lambda x: x ** 2
    else:
        func = lambda x: 1

    # 
    for i in range(1, n_points):
        x_prev = xi_mean[i - 1]        
        xi_mean[i] = x_prev + d * func(x_prev)
    
    xi_mean = xi_mean[:np.argmin(xi_mean < mz_end)]
    print("generating mean spectrum with {} sampling points...".format(len(xi_mean)))
    
    # "splat" warped peaks on sampling points to generate mean spectrum
    yi_mean = np.zeros(len(xi_mean))
    for i in range(0, len(spectra), stride):
        s_i = spectra[i]
        y_i = np.array(mx.splat_peaks(s_i, xi_mean, 3.0))
        yi_mean += y_i / tic[i]
    
    # use the 100 most intense mean spectrum peaks as references
    (xp, yp) = parabolic_centroid(xi_mean, yi_mean, 0.0)
    s_m = to_mx_peaks(xp, yp, sigma_1, len(spectra) + 100, instrument_type)
    
    print("generated mean spectrum")
    
    return s_m

def plot_range(fpath, mz0, mz1, ax, col, size=20, in_ppm = False):
    """Plot peak scatter within mz0 and mz1
    """
    ts = mx.get_triplets_range(fpath, mz0, mz1)
    
    if len(ts) == 0:
        return
    
    h = np.array([t.height for t in ts])
    t = np.array([t.index for t in ts])    
    mz = np.array([t.mz for t in ts])
    
    if in_ppm:
        x_m = np.mean(mz)
        x = (mz - x_m) / x_m * 1e6
        ax.scatter(x, t, size * h / np.max(h), c=col, alpha=0.5)   
    else:
        ax.scatter(mz, t, size * h / np.max(h), c=col, alpha=0.5)   


def triplets_to_ion_image(ts, meta):
    """Sum all intensities within mass bin for each pixel.
    returns: numpy matrix
    """
    s = np.zeros(len(meta))
    for t in ts:
        s[t.index] += t.height
        
    return to_ion_image(s, meta)


def plot_peak_matches(pms, exp_factor, pointsize, ax, **kwargs):
    """Plot
    """
    x = [pm[0].mz for pm in pms]
    y = [pm[1].mz - pm[0].mz for pm in pms]    
    h = np.array([np.power(pm[0].height * pm[1].height, exp_factor) for pm in pms])
    
    ax.scatter(x,y,h / np.max(h) * pointsize, **kwargs) 
    
    
def plot_warping(nodes, warping, ax):
    """Plot the recalibration function.
    """
    y = np.array([n.mz_shifts[w] for (n, w) in zip(nodes, warping)])
    x = np.array([n.mz for n in nodes])
    ax.plot(x, -y, '-o', c='tab:orange')