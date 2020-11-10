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




## Misc
