#! /usr/bin/env python3


"""

Will work with the integral of a rdf


"""

import numpy as np
import matplotlib.pyplot as plt
import pykbi


# read data from file

data = np.loadtxt("rdf1.txt")


rdf11 = pykbi.RDF(data[:,0],data[:,1], closed=True, name="rdf_11")
rdf22 = pykbi.RDF(data[:,0],data[:,2], closed=True, name="rdf_22")
rdf33 = pykbi.RDF(data[:,0],data[:,3], closed=True, name="rdf_33")


# we integrate them using the closed integration methods
rdf11.Integrate()
rdf22.Integrate()
rdf33.Integrate()

## and we read out the value of the integral. if no argument is given, this is done
# in the last point in the integral
rdf11.FindValues(position=(0.35,0.45))
rdf22.FindValues(position=(0.35,0.45))
rdf33.FindValues(position=(0.35,0.45))

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


rdf11.PlotRDF(ax1)
rdf22.PlotRDF(ax1)
rdf33.PlotRDF(ax1)

rdf11.PlotKBIInverse(ax2)
rdf22.PlotKBIInverse(ax2)
rdf33.PlotKBIInverse(ax2)

rdf11.PlotReadout(ax2)
rdf22.PlotReadout(ax2)
rdf33.PlotReadout(ax2)


#rdf.PlotKBI(ax2)
#rdf.PlotReadout(ax2)

#ax1.set_xlim(0.0, 5.0)

ax2.set_xlim(0.0, 1.0)

plt.show()
