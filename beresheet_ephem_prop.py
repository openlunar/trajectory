#!/usr/bin/env python3

"""
Script that takes one of the states from the BeresheetTrajectory.txt file and
propagates it in the ephemeris
"""
import numpy as np
from scipy.linalg import norm

from spiceypy import spiceypy as spice

from orbit import Orbit
from frames import rotate_z

import patched_conic as pc
from patched_conic import PatchedConic

from spice_loader import SpiceLoader
import propagate as prop
from propagate import Dynamics

from propagate.forces import gravity, j2_gravity

if __name__ == '__main__':

    loader = SpiceLoader()

    # Arbitrary arrival date at SOI (assuming model)
    arrival_date = '21JUN2022'
    arrival_time = spice.utc2et(arrival_date)

    x0 = np.array([-6.45306220e+06, -1.19390455e+06, -8.56868387e+04, 1.30926529e+03, -6.82322034e+03, -3.53025651e+03])
    moon0 = np.array([1.404200996426005e3])


    dynamics = Dynamics(fun_earth = gravity,
                        fun_moon  = gravity)
    t, x = prop.propagate_to(dynamics, init.depart_time,
                                          np.hstack((x0, )),
                                          PatchedConic.r_soi, arrival_time + 2 * 24 * 3600.0)

    print("{}: {}, {}".format(t, norm(x[0:3] - x[6:9]), PatchedConic.r_soi))
    #full_state = np.hstack((init.x_depart_post, init.x_moon_depart))
    #import pdb
    #pdb.set_trace()
    #x = prop.propagate_to(prop.dynamics, init.depart_time,
    #                      full_state,
    #                      arrival_time,
    #                      max_step = 600.0)
