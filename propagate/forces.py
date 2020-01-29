import numpy as np

from scipy.linalg import norm


def gravity(t, x, mu = None, label = None):
    f = np.zeros_like(x)
    r = norm(x[:3])
    
    f[0:3] = x[3:6]
    f[3:6] = x[0:3] * -mu / r**3
    return f


def j2_gravity(t, x, mu = None, j2 = None, r_m = None):
    f = np.zeros_like(x)
    r = norm(x[0:3])
    
    f[0:3] = x[3:6]
    pj = -1.5 * mu * j2 * r_m**2 / r**5
    
    f[3] = -mu * x[0] / r**3 + pj * x[0] * (1 - 5 * x[2]**2 / r**2)
    f[4] = -mu * x[1] / r**3 + pj * x[1] * (1 - 5 * x[2]**2 / r**2)
    f[5] = -mu * x[2] / r**3 + pj * x[2] * (3 - 5 * x[2]**2 / r**2)

    return f
