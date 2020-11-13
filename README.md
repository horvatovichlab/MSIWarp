# MSIWarp
**MSIWarp** is a flexible tool to perform mass alignment of Mass Spectrometry Imaging (MSI) spectra. A key feature of MSIWarp is its compatibility with centroid spectra.

## Installation
Clone the repository:
```
git clone --recursive https://github.com/horvatovichlab/MSIWarp.git
```
and build the project by typing
```
python3 setup.py install --user
```
in the root directory of the project.

CMake, a C++17 compliant compiler, and Python 3 must be installed to build MSIWarp. Furthermore, this project depends on the python packages Numpy, Matplotlib, and pyimzML (to interact with imzML files).

## Quick start
The following script aligns a list of spectra:
```python
import msiwarp as mx

spectra =  ... # code to load the unaligned spectra
reference_spectrum =  ... 

# setup the node placement parameters
params = mx.params_uniform(...)
epsilon = 1.0 # peak matching threshold, relative to peak width
n_cores = 4

# find an m/z recalibration function for each spectrum
recal_funcs = mx.find_optimal_warpings_uni(spectra, reference_spectrum, params, epsilon, n_cores)

# use the recalibration functions to warp the spectra
warped_spectra = [mx.warp_peaks_unique(s_i, r_i) for (s_i, r_i) in zip(spectra, recal_funcs)

# ... code to store the warped spectra

```

## Details
MSIWarp models a spectrum as a list of peaks ordered by *m/z*, where each peak has the following four attributes: 
1. spectrum index
2. *m/z*
3. height
4. sigma, the modelled width of the peak

### Initializing an MSIWarp spectrum
We can initialize an MSIWarp spectrum from a list of *m/z* and intensity values:
```python
import numpy as np
mzs = np.array([...]) # peak m/z values
hs = np.array([...]) # peak heights

import msiwarp as mx
s = [mx.peak(i, mz_i, h_i, 1.0) for i, (mz_i, h_i) in enumerate(zip(mzs, hs))]
```

### Setting up the warping nodes
A warping node has three parameters: an m/z position, a slack, and the number of evaluation points within the slack. The slack and the number of evaluation points determine the search space of the warping function and its resolution, respectively.

There are two fundamental options for warping node placement: (i) use the same warping nodes for all spectra, or (ii) place the warping nodes uniquely for each spectrum. If we choose (i), we must setup the warping nodes prior to searching for the optimal aligments. If we choose (ii), we must select a node placement method and set provide the corresponding parameters. 

The following code will initialize four warping nodes between 150 and 1050 *m/z*. 
```python
node_mzs = [150, 450, 750, 1050]
node_deltas = [0.015, 0.045, 0.075, 0.105] # slacks = node_deltas * n_steps
n_steps = 25 
nodes = mx.initialize_nodes(node_mzs, node_deltas, n_steps)
```

If we chose the second option, we can use either the density-based node placement function or the ... function.
```
# setup the parameters for the uniform node placement function
params = mx.params_uniform(mx.Instrument.Orbitrap, # each instrument type has its own relationship between peak width and m/z
                           n_steps, # same as above
                           n_peaks, # the number of peaks per segment
                           max_n_nodes, # maximum number of nodes
                           mz_begin, # 
                           mz_end, #
                           slack # the slack relative to peak width
                           )
```

## Example: align a centroid data set in the imzML format
Load a centroided imzML data set into RAM using [pyimzML](https://github.com/alexandrovteam/pyimzML):

```python
import msiwarp as mx

from pyimzml.ImzMLParser import ImzMLParser
from msiwarp.util.warp import to_mx_peaks

p = ImzMLParser(imzml_path)
spectra = []

for idx, coords in enumerate(p.coordinates):
    mzs, hs = p.getspectrum(idx)    
    spectra.append(to_mx_peaks(mzs, hs,
                               sigma_1, id = idx,
                               instrument_type = 'orbitrap'))

```

then warp the spectra:

```python
# choose a reference spectrum
i_r = 200
s_r = spectra[i_r]

print("warping spectra...")

import time
t0 = time.time()
optimal_moves = mx.find_optimal_spectra_warpings(spectra, s_r, nodes, epsilon)
t1 = time.time()
print("found optimal warpings in {:0.2f}s".format(t1 - t0))

t2 = time.time()
warped_spectra = [mx.warp_peaks(s_i, nodes, o_i) for (s_i, o_i) in zip(spectra, optimal_moves)]
t3 = time.time()
print("warped spectra in {:0.2f}s".format(t3 - t2))
```

    warping spectra...
    found optimal warpings in 77.72s
    warped spectra in 7.11s

and finally, store the warped spectra in a new imzML file:

```python
from pyimzml.ImzMLWriter import ImzMLWriter
from msiwarp.util.warp import to_mz, to_height

output_imzml = 'output.imzML'
with ImzMLWriter(output_imzml) as w:
    for s_i, coords in zip(warped_spectra, p.coordinates):
        # writes data to the .ibd file
        w.addSpectrum(to_mz(s_i), to_height(s_i), coords)
```

We can also store the data set in the *MSI triplet* format:
```python
fpath_triplets_raw = "..."
if mx.spectra_to_triplets(fpath_triplets_raw, spectra):
    print("wrote raw MSI triplets to file")

fpath_triplets_warped = "..."
if mx.spectra_to_triplets(fpath_triplets_warped, warped_spectra):
    print("wrote warped MSI triplets to file")
```
which enables fast queries of all data set peaks within a mass range. After generating our triplet files, we can easily plot mass scatters:

```python
from msiwarp.util.warp import plot_range

mz_ref = [281.249, 885.553, 886.556] # m/z locations of 
mass_tolerance = 3 # ppm

fig, ax = plt.subplots(1, 3, figsize=(12,4), sharey=True)

for i, mz_i in enumerate(mz_ref):
    d = mass_tolerance * mz_i / 1e6 # -+ mass_tolerance around reference mass 
    mz0 = mz_i - d
    mz1 = mz_i + d    
    
    plot_range(fpath_triplets_raw, mz0, mz1, ax[i], 'tab:cyan', 25, in_ppm=True)
    plot_range(fpath_triplets_warped, mz0, mz1, ax[i], 'tab:orange', 25, in_ppm=True)
    
    ax[i].set_facecolor((0.0, 0.0, 0.0))
    ax[i].set_title('m/z {:0.3f}'.format(mz_i))
    ax[i].set_xticks([-mass_tolerance, 0, mass_tolerance])
    ax[i].set_xlabel('relative shift (ppm)')
    
ax[0].set_ylabel('spectrum index')
```

![DESI MASS SCATTER](/docs/mass_scatter_desi.png)

## Misc
