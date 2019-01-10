#! /usr/bin/env python3


"""

Correct the rdf using the inverse N correction method.


"""

import numpy as np
import matplotlib.pyplot as plt
import pykbi


# read data from file

lt = 14.8245984505

data = np.loadtxt("rdf1.txt")

## load the first data
rdf11 = pykbi.RDF(data[:,0],data[:,1], closed=True, npart=1200, box_size=lt, eqint=True, name="rdf_11")
rdf22 = pykbi.RDF(data[:,0],data[:,2], closed=True, npart=600, box_size=lt, eqint=True, name="rdf_22")
rdf33 = pykbi.RDF(data[:,0],data[:,3], closed=True, npart=600, box_size=lt, eqint=True, name="rdf_33")
rdf12 = pykbi.RDF(data[:,0],data[:,4], closed=True, npart=1200, box_size=lt, eqint=False, name="rdf_12")
rdf13 = pykbi.RDF(data[:,0],data[:,5], closed=True, npart=1200, box_size=lt, eqint=False, name="rdf_13")
rdf23 = pykbi.RDF(data[:,0],data[:,6], closed=True, npart=600, box_size=lt, eqint=False, name="rdf_23")

## load the second set of data
data2 = np.loadtxt("rdf2.txt")

lt2=9.8313208864

rdf11_2 = pykbi.RDF(data2[:,0],data2[:,1], closed=True, npart=350, box_size=lt2, eqint=True, name="rdf_11_2")
rdf22_2 = pykbi.RDF(data2[:,0],data2[:,2], closed=True, npart=175, box_size=lt2, eqint=True, name="rdf_22_2")
rdf33_2 = pykbi.RDF(data2[:,0],data2[:,3], closed=True, npart=175, box_size=lt2, eqint=True, name="rdf_33_2")
rdf12_2 = pykbi.RDF(data2[:,0],data2[:,4], closed=True, npart=350, box_size=lt2, eqint=False, name="rdf_12_2")
rdf13_2 = pykbi.RDF(data2[:,0],data2[:,5], closed=True, npart=350, box_size=lt2, eqint=False, name="rdf_13_2")
rdf23_2 = pykbi.RDF(data2[:,0],data2[:,6], closed=True, npart=175, box_size=lt2, eqint=False, name="rdf_23_2")


## do the scaling
rdf11_invn = pykbi.CorrectInverseN(rdf11, rdf11_2)
rdf22_invn = pykbi.CorrectInverseN(rdf22, rdf22_2)
rdf33_invn = pykbi.CorrectInverseN(rdf33, rdf33_2)
rdf12_invn = pykbi.CorrectInverseN(rdf12, rdf12_2)
rdf13_invn = pykbi.CorrectInverseN(rdf13, rdf13_2)
rdf23_invn = pykbi.CorrectInverseN(rdf23, rdf23_2)


# we integrate them using the closed integration methods
rdf11_invn.Integrate()
rdf22_invn.Integrate()
rdf33_invn.Integrate()
rdf12_invn.Integrate()
rdf13_invn.Integrate()
rdf23_invn.Integrate()

## and we read out the value of the integral. if no argument is given, this is done
# in the last point in the integral
rdf11_invn.FindValues(position=(0.29,0.4))
rdf22_invn.FindValues(position=(0.29,0.4))
rdf33_invn.FindValues(position=(0.29,0.4))
rdf12_invn.FindValues(position=(0.29,0.4))
rdf13_invn.FindValues(position=(0.29,0.4))
rdf23_invn.FindValues(position=(0.29,0.4))

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


rdf11_invn.PlotRDF(ax1)
rdf22_invn.PlotRDF(ax1)
rdf33_invn.PlotRDF(ax1)
rdf12_invn.PlotRDF(ax1)
rdf13_invn.PlotRDF(ax1)
rdf23_invn.PlotRDF(ax1)


rdf11_invn.PlotKBIInverse(ax2)
rdf22_invn.PlotKBIInverse(ax2)
rdf33_invn.PlotKBIInverse(ax2)
rdf12_invn.PlotKBIInverse(ax2)
rdf13_invn.PlotKBIInverse(ax2)
rdf23_invn.PlotKBIInverse(ax2)

rdf11_invn.PlotReadout(ax2)
rdf22_invn.PlotReadout(ax2)
rdf33_invn.PlotReadout(ax2)
rdf12_invn.PlotReadout(ax2)
rdf13_invn.PlotReadout(ax2)
rdf23_invn.PlotReadout(ax2)


#ax1.set_xlim(0.0, 5.0)

ax2.set_xlim(0.0, 1.0)

plt.show()
