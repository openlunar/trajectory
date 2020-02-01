import numpy as np

from scipy.linalg import norm


def gravity(x, mu, gradient = False, **kwargs):
    """Returns a vector consisting of the time derivative of the state
    (the first six elements) and then the derivative of acceleration
    with respect to position (a 3x3 matrix) if requested (or None).

    """
    f = np.zeros(6)
    r = norm(x[:3])
    
    f[0:3] = x[3:6]
    f[3:6] = x[0:3] * -mu / r**3
    if gradient:
        return f, gravity_gradient(x[0:3], mu, r)
    else:
        return f, None


def j2_gravity(x, mu, gradient = False, j2 = None, r_eq = None):
    """Does the same thing as gravity() but also includes the J2
    perturbation (albeit not in the gravity gradient)."""
    f = np.zeros(6)
    r = norm(x[0:3])

    pj = 3 * mu * j2 * r_eq**2 / r**5
    l  = 5 * (x[2] / r)**2
    
    f[0:3] = x[3:6]
    f[3:6] = -mu * x[0:3] / r**3
    f[3:5] += pj * x[0:2] * (1 - l) * x[2]
    f[5]   += pj *   x[2] * (3 - l) * x[2]
    
    
    f[3] = -mu * x[0] / r**3 + pj * x[0] * (1 - 5 * x[0]**2 / r**2)
    f[4] = -mu * x[1] / r**3 + pj * x[1] * (1 - 5 * x[1]**2 / r**2)
    f[5] = -mu * x[2] / r**3 + pj * x[2] * (3 - 5 * x[2]**2 / r**2)

    if gradient:
        return f, gravity_gradient(x[0:3], mu, r)
    else:
        return f, None


def gravity_gradient(r, mu, rmag):
    rmag = norm(r[0:3])
    return np.outer(r,r) * (3.0 * mu / rmag**5) - np.identity(3) * (mu / rmag**3)
    
