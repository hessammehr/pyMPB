# pyMPB

This is a Python package for running bandstructure simulations of photonic crystals (PhCs) with MIT Photonic-Bands (MPB, http://ab-initio.mit.edu/wiki/index.php/MIT_Photonic_Bands). It creates a .ctl file, runs it, and does some postprocessing and plotting.

Features:
---------

* supports 2D and 3D simulations
* works with anisotropic dielectric materials
* supports MPB's cylinder and block objects (but easily extendable)
* rectangular and triangular lattices
* uses mpb-data (from MPB) and h5topng (from h5utils) to export PNG-Images of calculated fields
* exports simulation results (frequencies, group velocities) to .csv files
* plots band structure diagrams with matplotlib
* predefined simulations for 2D triangular lattice of holes and 3D triangular PhC slab
 

Dependencies:
-------------

* [MPB](http://ab-initio.mit.edu/wiki/index.php/MIT_Photonic_Bands)
* [H5utils](http://ab-initio.mit.edu/wiki/index.php/H5utils) (for png export of fields)
* [Python](https://www.python.org/)
* [numpy](https://pypi.python.org/pypi/numpy/)
* [matplotlib](http://matplotlib.org/) (for outputting diagrams)

 
Usage:
------
 
0. Add your material's refractive indices to data.py
0. Update defaults (programm calls etc.) in defaults.py
0. Start with predefined PhC simulations in phc_simulations.py
0. See examples in examples folder


"Disclaimer"
---------------

I am in a hurry - I need to finish my thesis. So my code is neither
nice nor clean. 
Meanwhile, I hope somebody finds this useful!

