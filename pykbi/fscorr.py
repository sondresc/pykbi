#! /usr/bin/env python3

"""
Perform finite size correction on the radial distribution fuction.

We have two methods:
    - InverseN- correction, as described by Kruger
    - van der Vegt correction, as described by ...

The InverseN correction method takes two RDF objects as input, and returns a new RDF object.

The van der Vegt correction takes onw RDF object, and returns a new RDF object.
"""

import numpy as np
import scipy
import pykbi.rdf as _rdf

__all__ = ["CorrectInverseN", "CorrectVanDerVegt"]


def CorrectInverseN(rdf1, rdf2):
    """
    Correct the rdf based on scaling between systems.
    """

    # check that the size of the bins are the same
    if abs(rdf1.r[1]-rdf2.r[1]) > 1e-4:
        print("CorrectInverseN correction requires rdf's with same resolution")
        return False

    # check that the rdfs are correctly defined
    if rdf1.npart is None or rdf2.npart is None:
        print("RDFs must have a defined number of particles")
        return False

    # check if the number of particles are different
    if rdf1.npart == rdf2.npart:
        print("Equal number of particles in each of the RDFs. ")
        print("Cannot do Inverse-N finite size correction when same number of atoms")
        return False

    # We do the inverse-N finte size correction
    if len(rdf1.r) > len(rdf2.r):
        rdfref = rdf2
        rdfext = rdf1
    else:
        rdfref = rdf1
        rdfext = rdf2

    bins = len(rdfref.r)

    correction_mask = (rdfref.gr - rdfext.gr[:bins]) / ((rdfext.npart/rdfref.npart) - 1.0)

    # build a new rdf-object
    out_rdf = _rdf.RDF(rdfref.r.copy(),
                       rdfref.gr-correction_mask,
                       npart=rdfref.npart,
                       box_size=rdfref.lt,
                       eqint=rdfref.eqint,
                       name=rdfref.name+" invN-corrected")

    return out_rdf



def CorrectVanDerVegt(rdf):
    """
    Using the van der Vegt correction method to correct the Finite size effect of
    a rdf.

    We only do this up to half the box size.

    """

    # check that number of particles is defined
    if rdf.npart is None:
        print("Cannot do VDV correction on this system. No particles defined")
        return False

    if rdf.lt is None or rdf.volume is None:
        print("Cannot do VDV correction on this system. ")
        print("RDF without box side lengths or volume")
        return False

    if rdf.eqint is None:
        print("Cannot do VDV correction on this system. ")
        print("RDF must have eqint value set.")
        return False


    # calculate the density of the component in the box. If it is a pair, the
    # first component should be used as the ref.
    rho_ref = rdf.npart/rdf.volume

    # kronecer delta, 1 if X-X, 0 otherwise. here we do an conversion from bool to int
    krondelta = int(rdf.eqint)

    c1 = rdf.npart * (1.0 - ((4.0 * np.pi * rdf.r**3 / 3.0) / rdf.volume))

    c2 = rho_ref * 4.0 * np.pi * scipy.integrate.cumtrapz(
        (rdf.gr[:] - 1.0) * rdf.r[:]**2,
        rdf.r[:])

    ## add the extra element at the beginning. We lose on array size when we integrate
    c2 = np.concatenate((np.array([0.0]), c2))

    c3 = (c1 / (c1 - c2 - krondelta))

    # build a new rdf-object
    out_rdf = _rdf.RDF(rdf.r.copy(),
                       rdf.gr*c3,
                       npart=rdf.npart,
                       box_size=rdf.lt,
                       eqint=rdf.eqint,
                       name=rdf.name)

    return out_rdf
