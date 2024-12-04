# Pyrf96

![Python3](https://img.shields.io/badge/python-3.x-brightgreen.svg)
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>

_Receiver function modelling_


This repository contains ctypes based Python wrappers for the program `RF.F90` from Shibutani et al. (1996) for forward modelling of seismic receiver functions of a plane wave impinging on a layered crustal Earth model.


## Installation
```
pip install git+https://github.com/inlab-geo/pyrf96
```

## Documentation

Essentially this is a single function, `rfcalc`. Here is the docstring:

```
Calculate synthetic Seismic receiver function waveforms for a given earth model, with optional addition of correlated noise (added in frequency domain). 

This is a slim Fortran wrapper around RF.F90 from Shibutani et al. (1996), which uses the Thomson-Haskell matrix formulation. Further details may be found.

Args:
    model (np.ndarray)               : Triplet defining layered model. meaning depends on mytpe, with shape (npts,3).
    sn (float,optional)              : Signal to noise ratio used to add correlated Gaussian noise to output.
    mtype (int, optional)            : Indicator for format of velocity model (default=0)
                                       model(1,i) is Vs velocity of layer i; 
                                       model(2,i) is vpvs ratio of layer i;
                                       mtype = 0 -> model(0,i) is the depth of Voronoi nuclei defining layer i;
                                       mtype = 1 -> model(0,i) is the thickness of layer i;
                                       mtype = 2 -> model(0,i) is depth of lower interface of layer i;
    fs (float, optional)             : Sampling frequency (default=25 samples/s)
    gauss_a (float, optional)        : Number 'a' defining the width of the gaussian filter in the deconvolution (default=2.5)
    water_c (float, optional)        : Water level used in deconvolution (default=0.0001)
    angle (float, optional)          : Angle in degrees of incoming teleseismic plane wave from vertical (default=35 degrees)
    time_shift (float, optional)     : Time shift before the first p pusle (default=5s)
    ndatar (int,optional)            : Number of time time steps of output waveform
    v60 (float,optional)             : P-wave velocity (km/s) needed to compute the ray parameter from angle (default=8.043 km/s)
    seed (int,optional)              : Random set for noise gneration

Returns:
    time (np.array, size ndatar )  : Time series time in seconds.
    wdata (np.array, size ndatar)  : The Receiver function amplitude.
```
## Example

```python
import numpy as np
import pyrf96 

vtype = 1 # Set up a velocity model in Voronoi cell format
velmod = np.zeros([7,3]) # nuclei can be unordered in depth
velmod[0] = [8.36, 3.25, 1.7]
velmod[1] = [17.2, 3.00, 1.7]
velmod[2] = [0.019, 2.51, 1.7]
velmod[3] = [19.7, 3.56, 1.7]
velmod[4] = [41.7, 4.23, 1.7]
velmod[5] = [14.3, 2.96, 1.7]
velmod[6] = [49.9, 4.59, 1.7]
    
time, amp = pyrf96.rfcalc(velmod,sn=0.25,mtype=vtype) # calculate receiver functions with correlated noise
```
A more detailed example of its usage is in `examples/RF_direct search demo.ipynb`

## Licensing
`pyrf96` is released as BSD-2-Clause licence


## Citations and Acknowledgments

> Shibutani, T., Kennett, B. and Sambridge, M.,  (1996) Genetic algorithm inversion for receiver functions with application to crust and uppermost mantle structure beneath Eastern Australia, Geophys. Res. Lett., 23 , No. 4, 1829-1832, 1996. 

Thanks to T. Shibutani (Disaster Prevention Research Institute, Kyoto Univ.) for creating the Fortran subroutine upon which this is based.