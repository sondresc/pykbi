# pykbi: Python module to work on radial distribution functions and Kirkwood-Buff integrals

Kirkwood-Buff theory bridges microscopic structure and macroscopic properties
of isotropic liquids. The classical version of the theory deals with open and infinite
systems\[[1](#KB1951)\], while the new correction allows us to work on closed and finite
systems\[[2](#Kruger2013)\]. The code contains modules to correct finite-size effects in Radial
Distribution functions (RDFs) from Molecular Simulations, and integrate RDFs
obtained from both open and closed ensembles. 

The methods in this module is included in the references below. `pykbi` works
with radial distribution functions obtained from various softwares, so it is up
to the user to read in and make radial distance and the pair correlation
function available as vectors. The integral can be done for both open and
closed systems.  In addition, it has included correction methods to correct for
finite-size effects.

This package accompanies the paper "Kirkwood-Buff integrals from molecular simulation"
\[[3](#dawass2019)\].

## Installation

The program can be installed from pip, using:

```bash
pip3 install pykbi
```

or it can be installed directly from github:

```bash
pip3 install git+git://github.com/sondresc/pykbi.git
```

or you can clone the code, and install manually:

```bash
git clone https://github.com/sondresc/pykbi
cd pykbi
python setup.py install
```

## References
1. <a name="KB1951" />[J. G. Kirkwood and F. P. Buff, *J. Chem. Phys.* **19**, 774(1951).](https://doi.org/10.1063/1.1748352)
1. <a name="Kruger2013" />[P. Kruger, S. K. Schnell, D. Bedeaux, S. Kjelstrup, T. J. H. Vlugt, J.-M. Simon, *J. Phys. Chem. Lett.* **4**, 2(2013).](https://doi.org/10.1021/jz301992u)
1. <a name="dawass2019" />[N. Dawass, P. Kruger, S. K. Schnell, J.-M. Simon, T. J. H. Vlugt, *Fluid Phase Equilibria* **486**, (2019).](https://doi.org/10.1016/j.fluid.2018.12.027)
