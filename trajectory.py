#!/usr/bin/env python3

"""This script should take a given trajectory optimization strategy /
model and, based on dates of interest, should produce a starting state
(inertial position and velocity) for the trajectory.

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

from propagate.forces import gravity, j2_gravity, zero_gravity


class InitialState(object):
    def __init__(self, arrival_time):
        # Start with earth parking orbit
        leo = Orbit.circular(PatchedConic.mu_earth, 6378136.6 + 185000.0)

        # Get the state of the moon at arrival so we know how far out
        # we need to go and how fast.
        x_moon_arrive = spice.spkez(301, arrival_time, 'J2000', 'NONE', 399)[0] * 1000.0
        
        # Produce patched conic
        x, pcx = pc.optimize_deltav(np.array([49.9 * np.pi/180.0,
                                              leo.v + 3200.0]),
                                    1837400.0, leo.r, leo.phi,
                                    norm(x_moon_arrive[0:3]),
                                    norm(x_moon_arrive[3:6]),
                                    conjugate = True)

        depart_time = arrival_time - pcx.tof
        free_flight_sweep_angle = pcx.nu1 - pcx.nu0
        
        # Get state of moon at departure so we can figure out the
        # plane of our trajectory.
        x_moon_depart = spice.spkez(301, depart_time,  'J2000', 'NONE', 399)[0] * 1000.0

        # Get earth--moon frame at SOI arrival time
        rm0 = x_moon_depart[:3]
        rm0hat = rm0 / norm(rm0)
        rm1 = x_moon_arrive[:3]
        rm1hat = rm1 / norm(rm1)
        hhat = np.cross(rm0hat, rm1hat)
        hhat /= norm(hhat)
        T_eci_to_pqw = spice.twovec(rm1hat, 1, hhat, 3)
        
        # Get directions to initial and final vectors
        r1hat = T_eci_to_pqw.T.dot(rotate_z(pcx.gam1).dot(np.array([1.0, 0, 0])))
        r0hat = T_eci_to_pqw.T.dot(rotate_z(pcx.gam1 - free_flight_sweep_angle).dot(np.array([1.0, 0.0, 0.0])))
        v0hat = np.cross(hhat, r0hat)

        # post delta-v state:
        r0 = r0hat * leo.r
        v0 = v0hat * x[1]
        # pre delta-v state:
        v0m = v0hat * leo.v

        # arrival state:
        # r1 = r1hat * pcx.arrive.r
        # v1 = ?

        self.free_flight_sweep_angle = free_flight_sweep_angle
        self.depart_time = depart_time
        self.arrival_time = arrival_time
        self.x_depart_post = np.hstack((r0, v0))
        self.x_depart_pre  = np.hstack((r0, v0m))
        self.x_moon_depart = x_moon_depart
        self.x_moon_arrive = x_moon_arrive
        self.deltav        = v0 - v0m
        #self.r1            = pcx.arrive.r
        #self.x_arrive = np.hstack((r1, v1))

    @property
    def r0(self):
        return self.x_depart_pre[:3]

    @property
    def v0_post(self):
        return self.x_depart_post[3:6]

    @property
    def v0_pre(self):
        return self.x_depart_pre[3:6]

    
def disperse(deltav,
             sigma_dtheta = 0.0,
             sfe          = 0.0):
    import pyquat as pq
    import pyquat.random as pqr
    #axis = axis_generator(**axis_generator_kwargs)

    # STUB: FIXME'


if __name__ == '__main__':

    loader = SpiceLoader()

    # Arbitrary arrival date at SOI (assuming model)
    arrival_date = '21JUN2022'
    arrival_time = spice.utc2et(arrival_date)
    

    init = InitialState(arrival_time)
    print("Initial ECI state is {}".format(init.x_depart_pre))


    # Can change gravity to j2_gravity, or alternatively to
    # zero_gravity if you want it to be off for a body.
    dynamics = Dynamics(fun_earth = gravity,
                        fun_moon  = gravity)
    x0 = np.hstack((init.x_depart_post, init.x_moon_depart, np.identity(6).reshape(36)))

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    axes = fig.add_subplot(111, projection='3d')
    
    t, x, Phi = prop.propagate_to_periselene(dynamics, init.depart_time, x0,
                                             arrival_time + 2 * 24 * 3600.0,
                                             plot = True,
                                             axes = axes)
    r_moon = x[0:3] - x[6:9]
    
    x01 = np.array(x0)
    x01[0] += 1.0
    
    t1, x1, Phi1 = prop.propagate_to(dynamics, init.depart_time, x01, t,
                                     plot = True,
                                     axes = axes)
    
    print("{}: {}, {}".format(t, norm(x[0:3] - x[6:9]), PatchedConic.r_soi))


    plt.show()

    import pdb
    pdb.set_trace()
