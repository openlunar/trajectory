#!/usr/bin/env python3

"""References:

* Arthur Gagg Filho, L., & da Silva Fernandes, S. 2016. Optimal round
  trip lunar missions based on the patched-conic approximation.
  Computational and Applied Mathematics, 35(3),
  753–787. https://doi.org/10.1007/s40314-015-0247-y

"""

import numpy as np

np.seterr(divide='raise', invalid='raise')

#from scipy.optimize import newton
from scipy.linalg import norm

from orbit import Orbit

def rotate_2d(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])

class PatchedConic(Orbit):
    # Physical constants of earth--moon system
    OMEGA = 2.649e-6 # r/s (mean)
    
    R = 1737400.0 # m -- radius of moon
    mu = 4.902800066163796e12 # m^3/s^2 moon gravitational constant
    mu_earth = 3.986004354360959e14 # m^3/s^2 earth gravitational constant
    r_soi = (mu / mu_earth)**0.4 * 384402000.0 # use mean distance for SOI calculation
    
    def __init__(self, depart, arrive,
                 lam1      = 0.0,
                 rf        = 1837000.0,
                 D         = 384402000.0,
                 V         = 2.649e-6 * 384402000.0):
        """Construct a planar patched conic approximation of an earth--moon
        transfer.

        Args:
          depart  earth-centered orbit at time of departure (post delta-v;
                  delta-v is assumed to be relative to a circular orbit of
                  the same radius)
          arrive  earth-centered orbit at SOI intercept
          lam1    spacecraft phase angle at arrival (angle between vector to
                  sphere-of-influence intercept and moon--earth
                  vector)
          rf      final desired radius of lunar orbit
          D       lunar distance at SOI arrival (defaults to the mean)
          V       lunar velocity at SOI arrival (defaults to the mean)
        """
        self.D = D
        self.V = V

        self.depart  = depart
        self.arrive  = arrive
        self.lam1    = lam1
        self.rf      = rf

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
        if arrive.a <= 0:
            raise ValueError("expected elliptical trajectory")
        tof  = np.sqrt(arrive.a**3 / arrive.mu) * ((E1 - e * sE1) - (E0 - e * sE0))


        # Get true anomalies at departure and arrival; we'll need
        # these to get departure phase angle.
        cnu0 = depart.cos_nu
        cnu1 = arrive.cos_nu
        nu0  = np.arccos(cnu0)
        nu1  = np.arccos(cnu1)

        # Get phase angle at arrival
        sg1 = np.clip((self.r_soi / arrive.r) * np.sin(lam1), -1.0, 1.0) # Eq. 4
        gam1 = np.arcsin(sg1) # gam1 is opposite lam1 in the arrival
                              # triangle (see Fig 1); phase angle at
                              # arrival
        
        # gam0 is the phase angle at departure
        # Note: This is approximate, since it depends on a fixed OMEGA
        gam0 = nu1 - nu0 - gam1 - self.OMEGA * tof

        # Eq. 9: velocity relative to moon when we reach SOI
        phi1 = arrive.phi # computed by Orbit.at()
        v1 = arrive.v
        v2 = np.sqrt(v1**2 + V**2 - 2.0 * v1 * V * np.cos(phi1 - gam1))

        # Compute miss angle of hyperbolic trajectory
        seps2 = np.clip( (V * np.cos(lam1) - v1 * np.cos(lam1 + gam1 - phi1)) / -v2, -1.0, 1.0 )
        eps2  = np.arcsin(seps2)

        # Eq. 10: Get selenocentric flight path angle 
        # right-hand side:
        tan_lam1_pm_phi2 = - v1 * np.sin(phi1 - gam1) / (V - v1 * np.cos(phi1 - gam1))
        phi2 = np.arctan(tan_lam1_pm_phi2) - lam1 # flight path angle

        # Parameters for Orbit class
        self.r   = self.r_soi
        self.v   = v2
        self.phi = phi2
        self.eps = eps2

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
        self.af   = self.r / (2.0 - Q2)

        # Position and velocity at perilune
        self.rpl  = self.af * (1.0 - self.ef)
        self.vpl  = np.sqrt((self.mu * (1.0 + self.ef)) / (self.af * (1.0 - self.ef)))

        self.deltav1 = np.abs(self.depart.v - np.sqrt(self.mu_earth / self.depart.r))
        self.deltav2 = self.vpl - self.vf
        self.tof     = tof

        # Objectives
        self.g = self.rf - self.rpl
        self.f = self.deltav1 + self.deltav2
        self.P = self.g**2

    def plot(self, alpha = 1.0, ax = None, v_scale = 100000.0):

        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        from scipy.linalg import norm

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

            # Plot moon, earth
            moon = Circle( (x.D, 0.0), 1737400.0, fc='grey', ec='grey', alpha=0.5)
            earth = Circle( (0.0, 0.0), 6378136.6, fc='blue', ec='blue', alpha=0.5)
            ax.add_patch(moon)
            ax.add_patch(earth)
            
            # Plot orbit of moon
            moon_orbit = Circle( (0.0, 0.0), self.D, fill=False, fc=None, alpha=0.5)
            ax.add_patch(moon_orbit)

            # Plot lunar SOI
            soi = Circle( (self.D, 0.0), self.r_soi, fill=False, fc=None, alpha=0.5)
            ax.add_patch(soi)

        # Plot from earth to intercept point, intercept point to moon
        # Find intercept point, call it r1
        r1 = rotate_2d(self.gam1).dot(np.array([self.arrive.r, 0.0]))
        ax.plot([0.0, r1[0], self.D],
                [0.0, r1[1], 0.0], c='k', alpha=alpha)

        # Plot moon-relative velocity
        r2_moon = rotate_2d(-self.lam1).dot(np.array([-self.r_soi, 0.0]))
        r2_earth = r2_moon + np.array([x.D, 0.0])
        v2_earth = rotate_2d(self.eps).dot(r2_moon / norm(r2_moon)) * -self.v * v_scale
        ax.plot([r2_earth[0], r2_earth[0] + v2_earth[0]],
                [r2_earth[1], r2_earth[1] + v2_earth[1]], c='r', alpha=alpha)

        # Plot earth-relative velocity
        v1 = rotate_2d(np.pi/2 - self.arrive.phi).dot(r1 / norm(r1)) * self.arrive.v * v_scale
        ax.plot([r1[0], r1[0] + v1[0]], [r1[1], r1[1] + v1[1]], c='g', alpha=alpha)

        vm = np.array([0.0, -self.V]) * v_scale
        ax.plot([r1[0] + v1[0], r1[0] + v1[0] + vm[0]], [r1[1] + v1[1], r1[1] + v1[1] + vm[1]], c='b', alpha=alpha)
        
        return ax


class PatchedConicGradients(object):
    def __init__(self, patched_conic):
        
        # Setup some shorthand notations
        v0 = patched_conic.depart.v
        v1 = patched_conic.arrive.v
        v2 = patched_conic.v
        vM = patched_conic.V
        
        Q2 = patched_conic.Q

        phi0 = patched_conic.depart.phi
        phi1 = patched_conic.arrive.phi
        gam1 = patched_conic.gam1

        phi2 = patched_conic.phi
        cphi2 = np.cos(phi2)
        sphi2 = np.sin(phi2)

        cphi1 = np.cos(phi1)
        sphi1 = np.sin(phi1)
        tphi1 = np.tan(phi1)
        cgam1 = np.cos(gam1)
        cpmg1 = np.cos(phi1 - gam1)
        spmg1 = np.sin(phi1 - gam1)

        ef = patched_conic.ef
        af = patched_conic.af

        lam1  = patched_conic.lam1

        mu      = patched_conic.depart.mu
        mu_moon = patched_conic.mu
        r0    = patched_conic.depart.r
        r1    = patched_conic.arrive.r
        r2    = patched_conic.r
        D     = patched_conic.D
        slam1 = np.sin(lam1)
        clam1 = np.cos(lam1)
        h     = patched_conic.arrive.h

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
        self.daf_dv0     = (2 * af**2 * v2 * self.dv2_dv0) / mu_moon
        self.drpl_daf    = 1.0 - ef # rpl = r_perilune
        self.drpl_def    = af
        self.drpl_dv0    = self.drpl_daf * self.daf_dv0 - self.drpl_def * self.def_dv0

        # Optimization state for Newton's method / restoration
        self.x = np.array([lam1, v0])



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

        self.df_dx          = np.array([[self.df_dlam1, self.df_dv0]]).T
        self.dg_dx          = np.array([[self.dg_dlam1, self.dg_dv0]]).T


        # SGRA computations
        dgdx = self.dg_dx # phi_x
        dfdx = self.df_dx # f_x
        P = 1.0 / dgdx.T.dot(dgdx)

        self.lam = (-dgdx.T.dot(dfdx) / (dgdx.T.dot(dgdx)))[0]
        
        # Compute derivative of merit function F
        self.dF_dx = dfdx + dgdx.dot(self.lam)

        # Compute stopping condition
        self.Q = self.dF_dx.T.dot(self.dF_dx)[0,0]

            


    #def dF_dx(self, lam):
    #    """Augmented function F derivative (for SGRA)"""
    #    return self.df_dx + self.dg_dx.dot(lam)

    #def F(self, lam):
    #    """Augmented function for SCGRA."""
    #    return self.f + self.g * lam

    #def Q(self, lam):
    #    """Stopping condition for SCGRA"""
    #    Fx = self.dF_dx(lam)
    #    return Fx.T.dot(Fx)[0,0]


def init_patched_conic(x,
                       rf   = 1837400.0,
                       r0   = 6371000.0 + 185000.0,
                       phi0 = 0.0,
                       D    = 384402000.0,
                       V    = 2.649e-6 * 384402000.0):
    """Generic objective function. x is [lam1, v0]."""
    lam1      = x[0]
    v0        = x[1]
    r_soi     = PatchedConic.r_soi
    r1        = np.sqrt(D**2 + r_soi**2 - 2.0 * D * r_soi * np.cos(lam1))
    depart    = Orbit(PatchedConic.mu_earth, r0, v0, phi0)
    intercept = depart.at(r1, sign='+')
    if np.isnan(intercept.v):
        raise ValueError("expected radius is not reached")
    elif depart.energy >= 0:
        raise ValueError("expected elliptical orbit")

    return PatchedConic(depart, intercept, lam1 = lam1, rf = rf)

def patched_conic_g(x, *args):
    """Objective function for achieving the final radius (also the constraint)."""
    return init_patched_conic(x, *args).g


def patched_conic_nr_g(x, lam1, *args):
    """Input for newton-raphson iteration"""
    return patched_conic_g(np.array([lam1, x]), *args)

def patched_conic_dg_dv0(x, lam1, *args):
    pcx = init_patched_conic(np.array([lam1, x]), *args)
    dx = PatchedConicGradients(pcx)
    print("g      = {}".format(pcx.g))
    print("dg_dv0 = {}".format(dx.dg_dv0))
    return dx.dg_dv0

def patched_conic_g_dg_dv0(x, lam1, *args):
    pcx = init_patched_conic(np.array([lam1, x]), *args)
    dx = PatchedConicGradients(pcx)
    print("g      = {}".format(pcx.g))
    print("dg_dv0 = {}".format(dx.dg_dv0))
    return pcx.g, dx.dg_dv0


def patched_conic_f(x, *args):
    """Objective function for minimizing total delta-v."""
    return init_patched_conic(x, *args).f
    

def patched_conic_df_dg(x, *args):
    pc = init_patched_conic(x, *args)
    dpc = PatchedConicGradients(pc)
    return (dpc.df_dv0, dpc.dg_dv0)

def dPsi_dalpha(alpha, x, *args):
    lam = args[-2]
    p = args[-1].reshape(2)
    args = args[0:-2]
    try:
        pcy = init_patched_conic(x - p * alpha, *args)
        dpcy = PatchedConicGradients(pcy)
        return -dpcy.Q
    except (ValueError, FloatingPointError):
        return np.float('nan')

    
def Psi(alpha, x, *args):
    lam = args[-2]
    p = args[-1].reshape(2)
    args = args[0:-2]
    try:
        pcy = init_patched_conic(x - p * alpha, *args)
        #dpcy = PatchedConicGradients(pcy)
        #lam = dpcy.lam
        
        #print("alpha = {}\tlam = {}\tf = {}\tg = {}".format(alpha, lam, pcy.f, pcy.g))
        #print("{}: {}".format(alpha, pcy.f + pcy.g * lam))
        return pcy.f + pcy.g * lam
    except (ValueError, FloatingPointError):
        return np.float('nan')


def Psi_dPsi_dalpha(alpha, x, *args):
    lam = args[-2]
    p = args[-1].reshape(2)
    args = args[0:-2]
    try:
        pcy = init_patched_conic(x - p * alpha, *args)
        Psi = pcy.f + pcy.g * lam
        dPsi_dalpha = dpcy.df_dx + dpcy.dg_dx * lam
        return (Psi, dPsi_dalpha)
        
    except (ValueError, FloatingPointError):
        #print("{}: nan".format(alpha))
        return (np.float('nan'), np.float('nan'))


def find_gradient(x, *args, conjugate = False,
                  dfdx         = None,
                  dgdx         = None,
                  dFdx2_prev   = None,
                  p_prev       = None,
                  alphatol     = 1e-8,
                  alphabracket = [1e-11, 0.1],
                  maxiter      = 100,
                  plot         = True):

    from scipy.optimize import minimize_scalar
    
    # Compute penalty parameter
    lam = (-dgdx.T.dot(dfdx) / (dgdx.T.dot(dgdx)))[0,0]

    # Compute gradient of merit function. Use this as the direction
    # for the alpha search.
    dFdx = (dfdx + dgdx * lam).flatten()
    if conjugate:
        dFdx2 = dFdx.T.dot(dFdx)
        if p_prev is None:
            p_prev = np.zeros_like(x)
            gamma  = 0.0
        else:
            gamma = dFdx2 / dFdx2_prev
        p = dFdx + p_prev * gamma
    else:
        p = dFdx
        dFdx2 = dFdx.T.dot(dFdx)

    # Find optimal alpha, where dPsi/dalpha is 0
    #
    # Use 'golden', because this appears to be the only method that
    # tolerates nans.
    alpha = newton(Psi,
                    alphabracket[0], (x, *args, lam, p),
                    tol     = alphatol,
                    maxiter = maxiter,
                    disp    = True,
                    minimize = True)

    if plot:
        try:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            da = (alpha - alphabracket[0]) / 50.0
            aa = np.arange(0.0, alpha + da, da)
            dpsi = np.zeros_like(aa)
            psi = np.zeros_like(dpsi)
            for ii in range(0, len(aa)):
                dpsi[ii] = dPsi_dalpha(aa[ii], x, *args, lam, p)
                psi[ii]  = Psi(aa[ii], x, *args, lam, p)
            ax.scatter(aa, psi)
            #ax.scatter(aa, dpsi)
            ax.axvline(alpha, c='r')
            plt.show()
        except ZeroDivisionError:
            pass
    


    return alpha, p, dFdx2

def find_restore_step(y, *args, maxiter=100, disp=True):
    pcy = init_patched_conic(y, *args)
    xt = y + 0.0
    pcxt = pcy
    dpcxt = PatchedConicGradients(pcxt)

    k = 1.0
    P_xt = float('inf')

    jj = 0
    while P_xt > pcy.P:
        #print("find_restore_step: {}".format(k))

        dgdx = dpcxt.dg_dx.flatten()
        
        sigma = pcxt.g * (k / (dgdx.dot(dgdx)))
        dy = dgdx * -sigma

        y_test = xt + dy.flatten()
        if disp:
            print("{}: xt = {}, g = {}, dg/dx = {}, sigma = {}, dy = {}".format(jj, xt, pcxt.g, dgdx, sigma, dy))
        try:
            pcy_test  = init_patched_conic(y, *args)
            dpcy_test = PatchedConicGradients(pcy_test)
            xt        = y_test
            pcxt      = pcy_test
            dpcxt     = dpcy_test
            P_xt      = pcxt.P
        except FloatingPointError:
            #pcxt = init_patched_conic(xt, *args)
            #dpcxt = PatchedConicGradients(pcxt)
            pass
        
        k *= 0.5
        jj += 1

        if jj == maxiter:
            raise ValueError("exceeded max iterations (restore step search)")

    return xt, pcxt, dpcxt


def restoration(y, *args, tol=1e-5, maxiter=100, sigma_maxiter=100):
    xt, pcxt, dpcxt = find_restore_step(y, *args, maxiter = sigma_maxiter)

    try:
        YS
    except NameError:
        YS = []
    
    ii = 0
    while pcxt.P > tol:
        y = xt

        YS.append(y)
            
        xt, pcxt, dpcxt = find_restore_step(y, *args, maxiter = sigma_maxiter)

        if norm(xt - y) < tol:
            print("Skipping restoration (y == xt)")
            return xt, pcxt, dpcxt

        if ii >= maxiter:
            raise ValueError("exceeded max iterations (restoration phase)")
        ii += 1
        
    return xt, pcxt, dpcxt


def newton_eval(fun_fprime, x, *args, step = 1e-11):
    fdata = fun_fprime(x, *args)
    if type(fdata) == tuple:
        f, df_dx = fdata
    else:
        f = fdata
        f1 = fun_fprime(x - step, *args)
        f2 = fun_fprime(x + step, *args)
        df_dx = (f2 - f1) / (2 * step)
    return f, df_dx

def newton(fun_fprime, x, args,
           tol         = 5e-5,
           step        = 1e-6,
           maxiter     = 100,
           beta0       = 1.0,
           beta_factor = 0.5,
           disp        = False,
           minimize    = False):
    """Find a departure velocity which allows us to fulfill our perilune
    constraint to within gtol. Use Newton's method, but must be
    sensitive to overshooting resulting in nans.

    """
    beta = beta0
    
    f, df_dx = newton_eval(fun_fprime, x, *args)
    for ii in range(0, maxiter):
        
        if disp:
            print("{}: x = {}, f = {}, dfdx = {}".format(ii, x, f, df_dx))

        # Stop if we reach our desired constraint tolerance
        if not minimize and np.abs(f) <= tol: # newton's method
            break

        try:
            dx = -beta * (f / df_dx)
            print("dx = {}".format(dx))
        except FloatingPointError:
            if minimize:
                break
            else:
                raise ZeroDivisionError("zero derivative encountered")

        # Keep track of previous values
        fp = f
        df_dxp = df_dx
        xp = x

        # Get next value
        x  = xp + dx
        
        # If update fails due to nan, try again with a smaller beta.
        try:
            f, df_dx = newton_eval(fun_fprime, x, *args)
            if minimize:
                if np.abs(dx) < tol:
                    break
            if np.abs(f) < np.abs(fp):
                retry = False
            else:
                retry = True
        except FloatingPointError:
            retry = True
        except ValueError:
            retry = True

        if retry:
            f = fp
            df_dx = df_dxp
            x = xp
            beta *= beta_factor
        elif beta < beta0:
            beta /= beta_factor
            if beta > beta0:
                beta = beta0
        #    beta = beta0 # reset beta

    if ii >= maxiter-1:
        raise ValueError("exceeded max iterations")

    return x


def optimize_deltav(x, *args,
                    maxiter          = 100,
                    gtol             = 5e-5,
                    Ptol             = 1e-5,
                    Qtol             = 2e-15,
                    alphatol         = 1e-12,
                    conjugate        = False,
                    newton_maxiter   = 100,
                    alpha_maxiter    = 100,
                    gradient_maxiter = 100,
                    restore_maxiter  = 100,
                    sigma_maxiter    = 100,
                    plot_alpha       = False):

    # For this lambda1, find the v0 that meets our constraint (correct
    # perilune radius) using Newton's method. This will be our
    # starting point.
    root_v0 = newton(patched_conic_g_dg_dv0,
                     x[1],
                     (x[0], *args),
                     tol = gtol,
                     maxiter = newton_maxiter,
                     disp = False)
    x[1] = root_v0
    print("x0 = {}".format(x))

    try:
        XS.append(x)
        YS.append(x)
    except NameError:
        XS = []
        YS = []

    pcx = init_patched_conic(x, *args)
    dpcx = PatchedConicGradients(pcx)

    alpha = 0.03
    # Get step (alpha), direction (p), and merit function gradient (dFdx)
    alpha, p, dFdx2 = find_gradient(x, *args,
                                    dfdx      = dpcx.df_dx,
                                    dgdx      = dpcx.dg_dx,
                                    conjugate = conjugate,
                                    maxiter   = alpha_maxiter,
                                    alphabracket = [0.0, alpha],
                                    alphatol  = alphatol,
                                    plot      = plot_alpha)
    print("alpha = {}, p = {}".format(alpha, p))


    for ii in range(0, gradient_maxiter):
        # Gradient phase:
        dx = -alpha * p.flatten()
        print("dx = {}".format(dx))
        y  = x + dx
        print("y = {}".format(y))

        # Restoration phase
        try:
            xt, pcxt, dpcxt = restoration(y, *args, tol = Ptol, maxiter = restore_maxiter, sigma_maxiter = sigma_maxiter)
            fail = False
        except ValueError as e:
            fail = True

        if fail or pcxt.f >= pcx.f: # reduce stepsize and repeat
            print("Reducing alpha step size to {} ({})".format(alpha, fail))
            if not fail:
                print("\tQ = {}".format(dpcxt.Q))
            alpha *= 0.5

        else: # Proceed to next gradient phase (unless we're beneath Qtol)
            print("post-restoration: alpha = {}\tf = {}\tg = {}\tQ = {}".format(alpha, pcxt.f, pcxt.g, dpcxt.Q))

            pcx  = pcxt
            x    = xt
            dpcx = dpcxt

            XS.append(x)

            # Normally this should be the finishing condition, but it
            # ain't working for some reason. That's fine. We'll assume
            # we're done when alpha won't get any smaller. It seems to
            # work.
            if dpcx.Q <= Qtol: # If we get here, we're done
                return x, pcx
            
            print("\tQ = {}".format(dpcx.Q))


            alpha, p, dFdx2 = find_gradient(x, *args,
                                            dfdx      = dpcx.df_dx,
                                            dgdx      = dpcx.dg_dx,
                                            conjugate = conjugate,
                                            maxiter   = alpha_maxiter,
                                            alphatol  = alphatol,
                                            plot      = plot_alpha)
    

    print("Warning: exceeded max iterations (gradient phase)")
    
    return x, pcx

if __name__ == '__main__':

    D      = 384402000.0
    V      = 2.649e-6 * D
    leo    = Orbit.circular(PatchedConic.mu_earth, 6378136.6 + 185000.0) # earth parking

    XS = []
    YS = []


    optimize_deltav(np.array([49.9 * np.pi/180.0,
                              leo.v + 3200.0]),
                    1837400.0, leo.r, leo.phi, D, V,
                    conjugate = True)


    YS = np.vstack(YS)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("SGRA optimization progress")
    ax.set_xlabel("\lambda_1 (deg)")
    ax.set_ylabel("v_0 (m/s)")
    ax.plot(np.array(YS[:,0]) * 180/np.pi, np.array(YS[:,1]), alpha=0.7)
    ax.text(YS[0,0] * 180/np.pi, YS[0,1], s='initial')
    ax.text(YS[-1,0] * 180/np.pi, YS[-1,1], s='final')
    plt.show()
    
        
    #opt    = SGRA()
    #alpha  = 1.0
    #for x in opt.optimize_v0(x, verbose=True):
    #    ax = x.plot(alpha = alpha, ax = ax)
    #    alpha *= 0.9   
    #plt.show()

