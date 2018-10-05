#! /usr/bin/env python3


"""
This small program will generate a rdf based on the Oscillatory Decaying
Function, integrate this based on the normal Kirkwood-Buff equation, and read
out the value in various points along the lenght of the integral.

It will print the state of the system, and make several figures.

"""

import numpy as np
import matplotlib.pyplot as plt
import pykbi



# first generate a range of values, we use excessively many points
range_data = np.linspace(0.001, 250, 10000)

# then we generate the function, based on the input
gr = pykbi.odf(range_data, chi=2.0)

rdf = pykbi.RDF(range_data, gr, closed=False, name="ODF, chi=2.0")

# we integrate the function using the open integration method.
rdf.Integrate()

## and we read out the value of the integral. if no argument is given, this is done
# in the last point in the integral
rdf.FindValues()


# Then print the values
print("Value of KBI at end of integral for '{}': {}".format(rdf.name, rdf.ReturnKBI()))

# We then print the state of the system
rdf.PrintState()

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


rdf.PlotRDF(ax1)
rdf.PlotKBI(ax2)
rdf.PlotReadout(ax2)

ax1.set_xlim(0.0, 5.0)

# write the data to a json-file
rdf.SaveToJSON("ocf_open.json")

plt.show()
