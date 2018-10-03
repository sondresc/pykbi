#! /usr/bin/env python3


"""

Test the finite-size correction using van der Vegt correction


"""

import numpy as np
import matplotlib.pyplot as plt
import pykbi


# read data from file

lt = 14.8245984505

data = np.loadtxt("rdf1.txt")

rdf11 = pykbi.RDF(data[:,0],data[:,1], closed=False, npart=1200, box_size=lt, eqint=True, name="rdf_11")
rdf22 = pykbi.RDF(data[:,0],data[:,2], closed=False, npart=600, box_size=lt, eqint=True, name="rdf_22")
rdf33 = pykbi.RDF(data[:,0],data[:,3], closed=False, npart=600, box_size=lt, eqint=True, name="rdf_33")
rdf12 = pykbi.RDF(data[:,0],data[:,4], closed=False, npart=1200, box_size=lt, eqint=False, name="rdf_12")
rdf13 = pykbi.RDF(data[:,0],data[:,5], closed=False, npart=1200, box_size=lt, eqint=False, name="rdf_13")
rdf23 = pykbi.RDF(data[:,0],data[:,6], closed=False, npart=600, box_size=lt, eqint=False, name="rdf_23")

rdf11_vdw = pykbi.CorrectVanDerVegt(rdf11)
rdf22_vdw = pykbi.CorrectVanDerVegt(rdf22)
rdf33_vdw = pykbi.CorrectVanDerVegt(rdf33)
rdf12_vdw = pykbi.CorrectVanDerVegt(rdf12)
rdf13_vdw = pykbi.CorrectVanDerVegt(rdf13)
rdf23_vdw = pykbi.CorrectVanDerVegt(rdf23)


# we integrate them using the closed integration methods
rdf11_vdw.Integrate()
rdf22_vdw.Integrate()
rdf33_vdw.Integrate()
rdf12_vdw.Integrate()
rdf13_vdw.Integrate()
rdf23_vdw.Integrate()

## and we read out the value of the integral. if no argument is given, this is done
# in the last point in the integral
rdf11_vdw.FindValues(position=(0.2,0.3))
rdf22_vdw.FindValues(position=(0.2,0.3))
rdf33_vdw.FindValues(position=(0.2,0.3))
rdf12_vdw.FindValues(position=(0.2,0.3))
rdf13_vdw.FindValues(position=(0.2,0.3))
rdf23_vdw.FindValues(position=(0.2,0.3))

# Then print the values
#print("Value of KBI at end of integral for '{}': {}".format(rdf.name, rdf.ReturnKBI()))

# We then print the state of the system
#rdf.PrintState()

# Before we finally make two figures showing the RDF, the integral, and the
# points where the data has been extracted.
fig1 = plt.figure()
fig2 = plt.figure()

ax1 = fig1.add_subplot(1, 1, 1)
ax2 = fig2.add_subplot(1, 1, 1)

ax1.set_xlabel("r")
ax1.set_ylabel("g(r)")

ax2.set_xlabel("R")
ax2.set_ylabel("G(R)")


rdf11_vdw.PlotRDF(ax1)
rdf22_vdw.PlotRDF(ax1)
rdf33_vdw.PlotRDF(ax1)
rdf12_vdw.PlotRDF(ax1)
rdf13_vdw.PlotRDF(ax1)
rdf23_vdw.PlotRDF(ax1)


rdf11_vdw.PlotKBIInverse(ax2)
rdf22_vdw.PlotKBIInverse(ax2)
rdf33_vdw.PlotKBIInverse(ax2)
rdf12_vdw.PlotKBIInverse(ax2)
rdf13_vdw.PlotKBIInverse(ax2)
rdf23_vdw.PlotKBIInverse(ax2)

rdf11_vdw.PlotReadout(ax2)
rdf22_vdw.PlotReadout(ax2)
rdf33_vdw.PlotReadout(ax2)
rdf12_vdw.PlotReadout(ax2)
rdf13_vdw.PlotReadout(ax2)
rdf23_vdw.PlotReadout(ax2)


#ax1.set_xlim(0.0, 5.0)

ax2.set_xlim(0.0, 1.0)

plt.show()
