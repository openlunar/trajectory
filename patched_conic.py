#!/usr/bin/env python3

"""References:

* Arthur Gagg Filho, L., & da Silva Fernandes, S. 2016. Optimal round
  trip lunar missions based on the patched-conic approximation.
  Computational and Applied Mathematics, 35(3),
  753–787. https://doi.org/10.1007/s40314-015-0247-y

"""

import numpy as np

from orbit import Orbit

class PatchedConic(Orbit):
    # Physical constants of earth--moon system
    D = 384402000.0 # distance from earth to moon in m
    OMEGA = 2.649e-3 # r/s
    V = OMEGA * D # mean velocity of moon relative to earth in m/s
    
    R = 1737000.0 # m -- radius of moon
    mu = 4.9048695e12 # m^3/s^2 moon gravitational constant
    mu_earth = 3.986004418e14 # m^3/s^2 earth gravitational constant
    r_soi = (mu / mu_earth)**0.4 * D
    
    def __init__(self, depart, arrive,
                 lam1 = 0.0,
                 rf   = 1837000.0):
        """Construct a planar patched conic approximation of an earth--moon
        transfer.

        Args:
          depart  earth-centered orbit at time of departure (post delta-v)
          arrive  earth-centered orbit at SOI intercept
          lam1    spacecraft phase angle at arrival (angle between vector to
                  sphere-of-influence intercept and moon--earth
                  vector)
          rf      final desired radius of lunar orbit

        """

        self.depart = depart
        self.arrive = arrive
        self.lam1   = lam1
        self.rf     = rf

        # E: eccentric anomaly in the transit orbit; E0 is departure,
        #    E1 is SOI intercept
        cE0  = depart.cos_E
        cE1  = arrive.cos_E
        E0   = np.arccos(cE0)
        E1   = np.arccos(cE1)
        sE0  = np.sin(E0)
        sE1  = np.sin(E1)

        e    = arrive.e

        # Get duration from departure from earth (time 0) to SOI
        # intercept (time 1)
        tof  = np.sqrt(arrive.a**3 / arrive.mu) * ((E1 - e * sE1) - (E0 - e * sE0))


        # Get true anomalies at departure and arrival; we'll need
        # these to get departure phase angle.
        cnu0 = depart.cos_nu
        cnu1 = arrive.cos_nu
        nu0  = np.arccos(cnu0)
        nu1  = np.arccos(cnu1)

        # Get phase angle at arrival
        sg1 = np.clip((self.r_soi / arrive.r) * np.sin(lam1), -1.0, 1.0)
        gam1 = np.arcsin(sg1) # gam1 is opposite lam1 in the arrival
                              # triangle (see Fig 1); phase angle at
                              # arrival
        
        # gam0 is the phase angle at departure
        gam0 = nu1 - nu0 - gam1 - self.OMEGA * tof

        # Eq. 9: velocity relative to moon when we reach SOI
        v1 = arrive.v
        v2 = np.sqrt(v1**2 + self.V**2 - 2.0 * v1 * self.V * np.cos(arrive.phi - gam1))

        # Angle of selenocentric velocity relative to moon's center
        phi1 = arrive.phi

        # Compute miss angle of hyperbolic trajectory
        seps2 = np.clip( (self.V * np.cos(lam1) - v1 * np.cos(lam1 + gam1 - phi1)) / -v2, -1.0, 1.0 )
        eps2  = np.arcsin(seps2)

        # Eq. 10: Get selenocentric flight path angle 
        # right-hand side:
        tan_lam1_pm_phi2 = - v1 * np.sin(phi1 - gam1) / (self.V - v1 * np.cos(phi1 - gam1))
        phi2 = np.arctan(tan_lam1_pm_phi2) - lam1 # flight path angle

        # Parameters for Orbit class
        self.r   = self.r_soi
        self.v   = v2
        self.phi = phi2

        # Additional parameters
        self.Q   = self.r * self.v**2 / self.mu
        self.rf  = rf
        self.vf  = np.sqrt(self.mu / self.rf)
        self.gam0 = gam0
        self.E0   = E0
        self.nu0  = nu0
        self.t0   = 0.0
        self.gam1 = gam1
        self.E1   = E1
        self.nu1  = nu1
        self.t1   = tof

        # Calculate eccentricity and semimajor axis using Eqs 52--53.
        Q2 = self.Q
        self.ef   = np.sqrt(1.0 + Q2 * (Q2 - 2.0) * np.cos(phi2)**2)
        self.af   = orbit.r / (2.0 - Q2)

        # Position and velocity at perilune
        self.rpl  = self.af * (1.0 - self.ef)
        self.vpl  = np.sqrt((self.mu * (1.0 + self.ef)) / (self.af * (1.0 - self.ef)))

        self.compute_gradients()
        
    def compute_gradients(self):
        orbit = self
        
        # Setup some shorthand notations
        v0 = orbit.depart.v
        v1 = orbit.arrive.v
        v2 = orbit.v
        vM = orbit.V
        
        Q2 = orbit.Q

        phi0 = orbit.depart.phi
        phi1 = orbit.arrive.phi
        gam1 = orbit.gam1

        phi2 = orbit.phi
        cphi2 = np.cos(phi2)
        sphi2 = np.sin(phi2)

        cphi1 = np.cos(phi1)
        sphi1 = np.sin(phi1)
        tphi1 = np.tan(phi1)
        cgam1 = np.cos(gam1)
        cpmg1 = np.cos(phi1 - gam1)
        spmg1 = np.sin(phi1 - gam1)

        ef = self.ef
        af = self.af



         # Eq. 57--61
        self.dv1_dv0     = v0 / v1
        self.dphi1_dv0   = (v0 / v1 - v1 / v0) / (v1 * tphi1)
        self.dv2_dv0     = ((v1 - vM * cpmg1) / v2) * self.dv1_dv0 + ((v1 * vM * spmg1) / v2) * self.dphi1_dv0
        self.dphi2_dv1   = -vM * spmg1 / v2**2
        self.dphi2_dphi1 = (v1**2 - v1 * vM * cpmg1) / v2**2
        self.dphi2_dv0   = self.dphi2_dv1 * self.dv1_dv0 + self.dphi2_dphi1 * self.dphi1_dv0

        # Eq. 63--65
        self.def_dv2     = 2 * Q2 * (Q2 - 1.0) * cphi2**2 / (ef * v2)
        self.def_dphi2   = -Q2 * (Q2 - 2.0) * cphi2 * sphi2 / ef
        self.def_dv0     = self.def_dv2 * self.dv2_dv0 + self.def_dphi2 * self.dphi2_dv0
        self.daf_dv0     = (2 * af**2 * v2 * self.dv2_dv0) / self.mu
        self.drpl_daf    = 1.0 - ef # rpl = r_perilune
        self.drpl_def    = af
        self.drpl_dv0    = self.drpl_daf * self.daf_dv0 - self.drpl_def * self.def_dv0

        # Optimization state for Newton's method
        self.x = np.array([orbit.lam1, v0])

        # Additional shorthand variables needed for SGRA optimization
        mu      = orbit.depart.mu
        mu_moon = orbit.mu
        r0 = orbit.depart.r
        r1    = orbit.arrive.r
        r2    = orbit.r
        D     = orbit.D
        lam1  = orbit.lam1
        slam1 = np.sin(orbit.lam1)
        clam1 = np.cos(orbit.lam1)
        h     = orbit.arrive.h
        
        self.deltav1 = np.abs(v0 - np.sqrt(mu / r0))
        self.deltav2 = orbit.vpl - orbit.vf

        self.dv1_dlam1 = -mu * D * r2 * slam1 / (v1 * r1**3)    # Eq. 66
        self.dphi1_dlam1   = h * D * r2 * slam1 / (v1 * r1**3 * sphi1) - h * D * r2 * mu * slam1 / (v1**3 * r1**4 * sphi1) # Eq. 67
        self.dgam1_dlam1 = r2 * clam1 / (r1 * cgam1) - D * (r2 * slam1)**2 / (r1**3 * cgam1) # Eq. 68

        self.dv2_dlam1 = ((v1 - vM * cpmg1) * self.dv1_dlam1
                          + (v1 * vM * spmg1) * self.dphi1_dlam1
                          - (v1 * vM * spmg1) * self.dgam1_dlam1) / v2 # Eq. 69
        self.dphi2_dgam1 = (vM * v1 * cpmg1 - v1**2) / v2**2 # Eq. 71

        # Eq. 73: Note --- this is only one portion of this equation;
        # we aren't using the rest right now.
        self.dphi2_dlam1 = self.dphi2_dphi1 * self.dphi1_dlam1 + self.dphi2_dgam1 * self.dgam1_dlam1 + self.dphi2_dv1 * self.dv1_dlam1 - 1.0

        self.dQ2_dlam1 = 2 * r2 * v2 * self.dv2_dlam1 / mu_moon # Eq. 74
        self.daf_dQ2   = af / (2.0 - Q2) # Eq. 75
        self.def_dQ2   = (Q2 - 1.0) * cphi2**2 / ef # Eq. 76

        self.daf_dlam1 = self.daf_dQ2 * self.dQ2_dlam1 # Eq. 78
        self.def_dlam1 = self.def_dQ2 * self.dQ2_dlam1 + self.def_dphi2 * self.dphi2_dlam1 # Eq. 79
        self.drpl_dlam1 = (1.0 - ef) * self.daf_dlam1 - af * self.def_dlam1 # Eq. 80
        self.dg_drpl    = -1.0
        self.dg_dlam1   = self.dg_drpl * self.drpl_dlam1
        self.dg_dv0     = -self.drpl_dv0 # Eq. 88, but equation in paper should be negated
        self.ddeltav2_dvpl = 1.0
        self.dvpl_daf      = 0.5 * np.sqrt((mu_moon * (1.0 + ef)) / (af**3 * (1.0 - ef))) # Eq. 89 (negated)
        self.dvpl_def      = -np.sqrt(mu_moon / ((1.0 + ef) * af * (1.0 - ef)**3)) # Eq. 90 (negated)
        self.ddeltav2_def  = self.ddeltav2_dvpl * self.dvpl_def
        self.ddeltav2_daf  = self.ddeltav2_dvpl * self.dvpl_daf
        self.dvpl_dlam1     = self.dvpl_daf * self.daf_dlam1 + self.dvpl_def * self.def_dlam1
        self.ddeltav2_dlam1 = self.ddeltav2_dvpl * self.dvpl_dlam1
        self.df_dlam1 = self.ddeltav2_dlam1
        self.ddeltav1_dv0   = 1.0
        self.dvpl_dv0       = self.dvpl_def * self.def_dv0 + self.dvpl_daf * self.daf_dv0
        self.ddeltav2_dv0   = self.ddeltav2_dvpl * self.dvpl_dv0
        self.df_dv0         = self.ddeltav1_dv0 + self.ddeltav2_dv0
        self.f              = self.deltav1 + self.deltav2
        self.df_dx          = np.array([[self.df_dlam1, self.df_dv0]]).T
        self.dg_dx          = np.array([[self.dg_dlam1, self.dg_dv0]]).T
        self.g              = orbit.rf

        # SGRA computations
        dgdx = self.dg_dx # phi_x
        dfdx = self.df_dx # f_x
        P = 1.0 / dgdx.T.dot(dgdx)
        self.lam = -P.dot(dgdx.T.dot(dfdx))[0,0]
        self.P = self.g**2

    def dF_dx(self, lam):
        """Augmented function F derivative (for SGRA)"""
        return self.df_dx + self.dg_dx.dot(lam)

    def F(self, lam):
        """Augmented function for SCGRA."""
        return self.f + self.g * lam

    def Q(self, lam):
        """Stopping condition for SCGRA"""
        Fx = self.dF_dx(lam)
        return Fx.T.dot(Fx)[0,0]

class SGRAOptimizer(object):
    D = PatchedConic.D
    V = PatchedConic.V
    OMEGA = PatchedConic.OMEGA
    mu_moon = PatchedConic.mu
    mu_earth = PatchedConic.mu_earth
    r_soi = PatchedConic.r_soi
    
    def __init__(self, r0, v0, rf,
                 lam1           = 70.0 * np.pi/180.0,
                 dv0            = 3150.0,
                 phi0           = 0.0,
                 rpg            = 185000.0 + 6371400.0,
                 gtol           = 5e-8,
                 ftol           = 2e-15,
                 alphatol       = 1e-6,
                 et_soi         = None,
                 max_iterations = 100,
                 beta           = 1.0):
        self.rf             = rf
        self.lam1           = lam1
        self.dv0            = dv0
        self.r0             = r0
        self.phi0           = phi0
        self.gtol           = gtol
        self.ftol           = ftol
        self.beta0          = beta

        # Find the radius of encounter with the moon's SOI
        self.r1             = np.sqrt(self.D**2 + self.r_soi**2 - 2.0 * self.D * self.r_soi * np.cos(lam1))

        if not self.update(v0):
            raise ValueError("orbit is invalid, produces no constraint gradient")

    def update(self, v0):
        # Use v0 to find v1 with the departure trajectory.
        depart = Orbit(self.mu_earth, self.r0, v0, self.phi0)
        intercept = depart.at(self.r1, sign='+')

        # Now we can get a hyperbolic arrival orbit (moon-centered)
        x = PatchedConic(depart, intercept, lam1 = self.lam1, rf = self.rf)

        if np.isnan(x.dg_dv0):
            import pdb
            pdb.set_trace()
            return False
        else:
            self.v0 = v0
            self.lunar = x
            return True

    def optimize_v0(self, max_iterations = 100):
        x = self.lunar
        beta = self.beta0
        
        for ii in range(0, max_iterations):

            # Stop if we reach our desired constraint tolerance
            if np.abs(x.g) <= self.gtol:
                break

            dv0 = -beta * (x.g / x.dg_dv0)

            prev_v0 = self.v0
            v0 += dv0

            # If update fails, try again with a smaller beta. If it
            # passes, reset beta to the initial, and keep on going.
            if self.update(v0):
                x = self.lunar
                beta = self.beta0

                yield x
            else:
                beta *= 0.5
                v0 = prev_v0

    def optimize_deltav(self, max_iterations = 100):
        raise NotImplemented("FIXME")



if __name__ == '__main__':

    orbit  = Orbit.circular(PatchedConic.mu_earth, 185000.0 + 6371400.0) # earth parking
    opt    = SGRAOptimizer(orbit.r, orbit.v, 1937.0)
    for x in opt.optimize_v0():
        import pdb
        pdb.set_trace()
