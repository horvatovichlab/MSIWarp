# MSIWarp

## Installation
This project depends on Numpy, Matplotlib, pyimzml (to interact with imzML files)

To build the MSIWarp Python package, cmake, a C++17 compliant compiler, and python 3 (version?) must be installed

To build the project, open a terminal in the root folder and type
```
python3 setup.py install --user
```

## Aligning a data set

MSIWarp models a spectrum as a list of peaks, where each peak has the following four attributes: 
1. spectrum index
2. m/z
3. height
4. sigma, the modelled width of the peak

### Initializing spectra
To generate a MSIWarp spectrum from a list of m/z and intensity values:
```python
import numpy as np
mzs = np.array([...]) # peak m/z values
hs = np.array([...]) # peak m/z heights

import msiwarp as mx
s = [mx.peak(i, mz_i, h_i, 1.0) for i, (mz_i, h_i) in enumerate(zip(mzs, hs))]
```
To search for an optimal warping between a pair of spectra we 

### Setting up warping nodes
The slack is equal to the number of steps times the step size (node delta)
```
nodes = mx.initialize_nodes(node_mzs, node_deltas, n_steps)
```
### Interfacing with pyimzml (TODO: add link to pyimzml homepage)

Load a centroided imzml data set into RAM

```python
import msiwarp as mx

from pyimzml.ImzMLParser import ImzMLParser
from msiwarp.util.warp import to_mx_peaks

positions = []
spectra = []

p = ImzMLParser(imzml_path)
for idx, coords in enumerate(p.coordinates):
    positions.append(coords)
    
    # 
    mzs, hs = p.getspectrum(idx)    
    spectra.append(to_mx_peaks(mzs, hs,
                               sigma_1, id = idx,
                               instrument_type = 'orbitrap'))

```

Warp spectra

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

and store warped spectra in a new imzML file:

```python
# convert back to pyimzml format

from pyimzml.ImzMLWriter import ImzMLWriter
from msiwarp.util.warp import to_mz, to_height

output_imzml = fdir + 'output.imzML'

with ImzMLWriter(output_imzml) as w:
    for s_i, coords in zip(warped_spectra, positions):
        # writes data to the .ibd file
        w.addSpectrum(to_mz(s_i), to_height(s_i), coords)
```

![DESI MASS SCATTER](/docs/mass_scatter_desi.png)

## Misc
