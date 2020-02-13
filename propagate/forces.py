import numpy as np

from scipy.linalg import norm

def CR3BP(X, mu, gradient = False):
    Xdot = np.zeros(6)
    x,y,z,vx,vy,vz = X[:]

    mu1 = 1-mu
    r13 = ((x+mu)**2 + y**2 + z**2)**1.5 #distance from M1 to P cubed
    r23 = ((x-mu1)**2 + y**2 + z**2)**1.5 #distance from M2 to P cubed

    Xdot[0:3] = vx,vy,vz
    Xdot[3] = 2*vy + x - mu1*(x+mu)/r13 - mu*(x-mu1)/r23
    Xdot[4] = -2*vx + y - mu1*y/r13 - mu*y/r23
    Xdot[5] = -mu1*z/r13 - mu*z/r23
    if gradient:
        G = np.zeros((3,3))
        r15 = ((x+mu)**2 + y**2 + z**2)**2.5 #distance from M1 to P raised to the fifth power
        r25 = ((x-mu1)**2 + y**2 + z**2)**2.5 #distance from M2 to P raised to the fifth power

        G[0,0]      = 1 - mu1*(1/r13 - 3*(x+mu)**2/r15) - mu*(1/r23 - 3*(x-mu1)**2/r25);
        G[0,1]      = 3*mu1*(x+mu)*y/r15 + 3*mu*(x-mu1)*y/r25;
        G[0,2]      = 3*mu1*(x+mu)*z/r15 + 3*mu*(x-mu1)*z/r25;
        G[1,1]      = 1 - mu1*(1/r13 - 3*y**2/r15) - mu*(1/r23 - 3*y**2/r25);
        G[1,2]      = 3*mu1*y*z/r15 + 3*mu*y*z/r25;
        G[2,2]      = - mu1*(1/r13 - 3*z**2/r15) - mu*(1/r23 - 3*z**2/r25);

        G[1,0]      = G[0,1]; #Symmetric Matrix
        G[2,0]      = G[0,2];
        G[2,1]      = G[1,2];
        return Xdot, G
    else:
        return Xdot, None

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
