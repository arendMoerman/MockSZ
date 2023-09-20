import numpy as np

from MockSZ.Constants import Constants as ct

def getXYGrid(x, y):
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

def getXYZGrid(x, y, z):
    if isinstance(y, float) and isinstance(z, float) and not isinstance(x, float):
        X, Y, Z = np.mgrid[x[0]:x[-1]:x.size*1j, y:y:1j, z:z:1j]

    elif not isinstance(y, float) and isinstance(z, float) and isinstance(x, float):
        X, Y, Z = np.mgrid[x:x:1j, y[0]:y[-1]:y.size*1j, z:z:1j]

    elif isinstance(y, float) and isinstance(z, float) and isinstance(x, float):
        X = np.array([x])
        Y = np.array([y])
        Z = np.array([z])

    else:
        X, Y, Z = np.mgrid[x[0]:x[-1]:x.size*1j, y[0]:y[-1]:y.size*1j, z[0]:z[-1]:z.size*1j]
   
    return X, Y, Z
