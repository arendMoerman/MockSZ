"""!
@file
Utility functions for MockSZ.
"""

import numpy as np

def getXYGrid(x, y):
    """!
    Generate a grid of two (arrays of) values.

    @param x Float or array of size nx containing the value(s) along the abscissa.
    @param y Float or array of size ny containing the value(s) along the ordinate.
    
    @returns X C-style array of shape (nx, ny). If both x and y are float, returns a singleton.
    @returns Y C-style array of shape (nx, ny). If both x and y are float, returns a singleton.
    """

    if isinstance(y, float) and not isinstance(x, float):
        X, Y = np.mgrid[x[0]:x[-1]:x.size*1j, y:y:1j]

    elif not isinstance(y, float) and isinstance(x, float):
        X, Y = np.mgrid[x:x:1j, y[0]:y[-1]:y.size*1j]

    elif isinstance(y, float) and isinstance(x, float):
        X = np.array([x])
        Y = np.array([y])

    else:
        X, Y = np.mgrid[x[0]:x[-1]:x.size*1j, y[0]:y[-1]:y.size*1j]
   
    return X, Y

