
"""
Oscillatory decaying function.

This is the same function as is used in Kruger et al. J. Phys. Chem. Lett.
2013, 4, 235-238 (dx.doi.org/10.1021/jz301992u)

The function acts as a radial distribution function, being zero initally, and
oscilating around 1.  The value r is the radial distance, chi is the intensity
of the fluctuations, and sigma is the width of the fluctuations.

"""

import numpy as _np

__all__ = ["odf"]

def odf(radius, chi, sigma=1.0):
    """
    Calculate the fluctuating decaying function

    """
    cut = 19.0 / 20.0

    rdata = radius/sigma

    function = lambda rin: 1.5/rin * \
            _np.exp((1.0 - rin)/chi) * _np.cos(2.0 * _np.pi * (rin - (21.0/20.0)))

    return _np.where((rdata >= cut), function(rdata), -1.0) + 1.0
