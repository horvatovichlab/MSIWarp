#!/usr/bin/env python
# coding: utf-8

import msiwarp as mx

from msiwarp.util.warp import to_mx_peaks, to_mz, to_height
from msiwarp.util.warp import generate_mean_spectrum, plot_range

import matplotlib.pyplot as plt
import numpy as np


# ----------------- setup -----------
fdir = 'datasets/orbitrap-desi/'

imzml_path = fdir + 'A67 CT S4-centroid.imzML'
fpath_triplets_raw = fdir + 'orbitrap_desi_triplets_raw.dat'
fpath_triplets_warped = fdir + 'orbitrap_desi_triplets_warped.dat'
fpath_dispersion_csv = fdir + 'results/dispersion_100.csv'
fpath_scatter = fdir + 'results/scatter'

# experiment settings
sigma_1 = 3.0e-7
epsilon = 1.0
instrument_type = 'orbitrap'


# --------- load spectra and meta data ----------
from pyimzml.ImzMLParser import ImzMLParser
from msiwarp.util.warp import to_mx_peaks

p = ImzMLParser(imzml_path)
spectra = []

for idx, coords in enumerate(p.coordinates):
    mzs, hs = p.getspectrum(idx)    
    spectra.append(to_mx_peaks(mzs, hs,
                               sigma_1, id = idx,
                               instrument_type = instrument_type))


tic = np.array([np.sum(to_height(s_i)) for s_i in spectra])
n_peaks = np.array([len(s_i) for s_i in spectra])


# ---------- use uniform node placement function ----------
n_steps = 33
n_peaks_min = 30
n_nodes = 8
mz_begin = 175
mz_end = 1000
slack = 2.0 * epsilon * sigma_1

params = mx.params_uniform(mx.Instrument.Orbitrap,
                           n_steps,
                           n_peaks_min,
                           n_nodes,
                           mz_begin,
                           mz_end,
                           slack)


# ---------- mean spectrum ----------
n_points = 2000000
s_m = generate_mean_spectrum(spectra, n_points, sigma_1,
                             mz_begin, mz_end, tic, instrument_type)

print("using mean spectrum as reference")
s_m_1000 = mx.peaks_top_n(s_m, 1000) # returns peak list sorted by intensity, not m/z
s_ref = sorted(s_m_1000, key=lambda peak: peak.mz)

s_m_100 = mx.peaks_top_n(s_m, 100)
mz_ref = np.sort(to_mz(s_m_100))


# ---------- warp spectra ----------
print("warping spectra...")

import time

n_cores = 8

t0 = time.time()
warping_funcs = mx.find_optimal_warpings_uni(spectra, s_ref, params, epsilon, n_cores)
t1 = time.time()
print("found optimal warpings in {:0.2f} seconds".format(t1 - t0))

t2 = time.time()
warped_spectra = [mx.warp_peaks_unique(s_i, r_i) for (s_i, r_i) in zip(spectra, warping_funcs)]
t3 = time.time()
print("warped spectra in {:0.2f}s".format(t3 - t2))


# ---------- save warped spectra as MSI triplets ----------
if mx.spectra_to_triplets(fpath_triplets_raw, spectra):
    print("wrote raw triplets to file")

if mx.spectra_to_triplets(fpath_triplets_warped, warped_spectra):
    print("wrote warped triplets to file")


# ---------- plot mass scatter around mean spectrum peaks ----------
mass_tolerance = 4 # ppm

import os

if not os.path.exists(fpath_scatter):    
    os.makedirs(fpath_scatter)
    print("made scatter plot output directory")
else:
    print("using existing scatter plot output directory")
    
fig, ax = plt.subplots(figsize=(6,6))
for i, mz_i in enumerate(mz_ref):
    d = mass_tolerance * mz_i / 1e6 # -+ 4 ppm around reference mass 
    mz0 = mz_i - d
    mz1 = mz_i + d    
    
    plot_range(fpath_triplets_raw, mz0, mz1, ax, 'tab:cyan', 25, in_ppm=True)
    plot_range(fpath_triplets_warped, mz0, mz1, ax, 'tab:orange', 25, in_ppm=True)
    
    ax.set_facecolor((0.0, 0.0, 0.0))
    ax.set_title('m/z {:0.3f}'.format(mz_i))
    ax.set_xticks([-mass_tolerance, 0, mass_tolerance])
    
    fig.savefig(fpath_scatter + '/mz_{}.png'.format(int(mz_i)), dpi=200)
    ax.cla()


# ---------- compute mass dispersions around mean spectrum ----------
from msiwarp.util.warp import dispersion_triplets
import pandas as pd

dispersion_raw = np.zeros(len(mz_ref))
dispersion_warped = np.zeros(len(mz_ref))
    
for i, mz_i in enumerate(mz_ref):
    d = mass_tolerance * mz_i / 1e6 # -+ 4 ppm around reference mass
    mz0 = mz_i - d
    mz1 = mz_i + d
    
    ts_raw = mx.get_triplets_range(fpath_triplets_raw, mz0, mz1)
    ts_warped = mx.get_triplets_range(fpath_triplets_warped, mz0, mz1)
    
    q = 0.0 # remove background signal
    if len(ts_raw) > 0:
        dispersion_raw[i] = dispersion_triplets(ts_raw,  q)
    if len(ts_warped) > 0:  
        dispersion_warped[i] = dispersion_triplets(ts_warped, q)


d = {'mz': mz_ref,
     'dispersion raw [ppm]': dispersion_raw,
     'dispersion warped [ppm]': dispersion_warped}

df = pd.DataFrame(d)
df.round(4).to_csv(fpath_dispersion_csv, index=False)

print('median mass dispersion raw: {:0.4f}'.format(np.median(dispersion_raw)))
print('median mass dispersion warped: {:0.4f}'.format(np.median(dispersion_warped)))