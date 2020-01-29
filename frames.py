import numpy as np
from scipy.linalg import norm

def rotate_x(t):
    """Rotation about the x axis by an angle t.

    Args:
        t  angle (radians)

    Returns:
        3x3 orthonormal rotation matrix.
    """
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, np.cos(t), -np.sin(t)],
                     [0.0, np.sin(t),  np.cos(t)]])


def rotate_y(t):
    """Rotation about the y axis by an angle t.

    Args:
        t  angle (radians)

    Returns:
        3x3 orthonormal rotation matrix.
    """
    return np.array([[np.cos(t), 0, np.sin(t)],
                     [0, 1, 0],
                     [-np.sin(t), 0, np.cos(t)]])

def rotate_z(t):
    """Rotation about the z axis by an angle t.

    Args:
        t  angle (radians)

    Returns:
        3x3 orthonormal rotation matrix.
    """
    return np.array([[np.cos(t), -np.sin(t), 0.0],
                     [np.sin(t),  np.cos(t), 0.0],
                     [0, 0, 1]])


def compute_T_inrtl_to_lvlh(x_inrtl):
    """Compute a state transformation from inertial frame to LVLH frame:

    * x is the local horizontal (v-bar)
    * y is the local out-of-plane
    * z is the local vertical (r-bar)

    The LVLH frame is a local frame, meaning that it is generally
    centered on your spacecraft (or other object of interest). It's
    useful for giving relative position or velocity of a nearby
    object, or for describing your covariance in terms that are useful
    for entry/descent/guidance (like, how well do I know my vertical
    velocity?).

    References:

    * Williams, J. (2014). LVLH Transformations.

    Args:
      x_inrtl  state vector (size 6 or 9: position and velocity are
               mandatory but acceleration is optional)

    Returns:
        Returns a 6x6 matrix for transforming a position and velocity
    state.

    """
    r = x_inrtl[0:3]
    v = x_inrtl[3:6]
    h = np.cross(r, v)

    # Normalize to get vectors for position, velocity, and orbital
    # angular velocity.
    r_mag = norm(r)
    v_mag = norm(v)
    
    r_hat = r / r_mag
    v_hat = v / v_mag
    
    h_mag = norm(h)
    h_hat = h / h_mag

    # Basis vector
    e_z  = -r_hat
    e_y  = -h_hat
    e_x  = np.cross(e_y, e_z) # Williams has a typo in the derivation here.
    e_x /= norm(e_x)

    de_z = (r_hat * np.dot(r_hat, v) - v) / r_mag

    # If acceleration is occurring, procedure is more complicated:
    if x_inrtl.shape[0] == 9:
        a     = x_inrtl[6:9]
        h_dot = np.cross(r, a)
        de_y  = (h_hat * np.dot(h_hat, h_dot) - h_dot) / h_mag
        de_x  = np.cross(de_y, e_z) + np.cross(e_y, de_z)
    else: # Simple two-body mechanics
        de_y = np.zeros(3)
        de_x = np.cross(de_z, h_hat)

    Til = np.vstack((e_x, e_y, e_z))
    dTil = np.vstack((de_x, de_y, de_z))

    T = np.vstack( (np.hstack( (Til,  np.zeros((3,3))) ),
                    np.hstack( (dTil, Til) )) )

    return T

    

def approx_T_pcpf_to_inrtl(t, omega = np.array([0.0, 0.0, 7.292115e-5])):
    return rotate_z(t * omega)
