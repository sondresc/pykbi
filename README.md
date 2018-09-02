# pykbi

## Python module to work on radial distribution functions and Kirkwood-Buff integrals

Kirkwood-Buff theory bridges microscopic structure and macroscopic properties
of isotropic liquids. The classical theory deals with open and infinite systems, while the 
new correction allows us to work on closed and finite systems. The code contains modules to correct finite-size 
effects in Radial Distribution functions (RDFs) from Molecular Simulations, and integrate RDFs obtained from 
both open and closed ensembles. 

References to the methods included in this program can be found in: 

:point_right: Krüger et al, ...


## Installation

The program can be installed from pip, using

```bash
pip install pykbi
```

or it can be installed directly from github:

```bash
pip install git+git://github.com/sondresc/pykbi.git
```

or you can clone the code, and install it directly on your computer:

```bash
git clone https://github.com/sondresc/pykbi
cd pykbi
python setup.py install
```
