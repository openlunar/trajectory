import numpy as np
from scipy.linalg import norm

class Orbit(object):

    @classmethod
    def elliptical(self, mu = 398600.436, rp = 3800.0, ra = 42000.0):
        a = (ra + rp) / 2.0
        eps = -mu / (2*a)
        vp = np.sqrt(2 * (eps + mu/rp) )

        return Orbit(mu, rp, vp, 0.0)

    @classmethod
    def circular(self, mu, r):
        return Orbit(mu, r, np.sqrt(mu / r))

    def __init__(self, mu, r, v, flight_path_angle = 0.0):
        self.mu = mu
        self.r = r
        self.v = v
        self.phi = flight_path_angle
    
    @property
    def energy(self):
        return self.v**2 / 2.0 - self.mu / self.r

    @property
    def h(self):
        return self.r * self.v * np.cos(self.phi)

    @property
    def a(self):
        return -self.mu / (2.0 * self.energy)

    @property
    def rp(self):
        return (1.0 - self.e) * self.a

    @property
    def ra(self):
        return (1.0 + self.e) * self.a

    @property
    def p(self):
        return self.h**2 / self.mu

    @property
    def vp(self):
        return self.v_at_radius(self.rp)

    @property
    def va(self):
        return self.v_at_radius(self.ra)

    @property
    def vinf(self):
        return np.sqrt(self.v**2 - 2 * self.mu / self.r)
        

    @property
    def e(self):
        return np.sqrt(np.clip(1.0 - self.p / self.a, 0.0, float('inf')))

    @property
    def period(self):
        return 2 * np.pi * np.sqrt(self.a**3 / self.mu)

    def at(self, r1, sign='+'):
        v1 = self.v_at_radius(r1)
        cos_phi1 = np.clip(self.h / (r1 * v1), -1.0, 1.0)

        if sign == '+':
            phi1 = np.arccos(cos_phi1)
        else:
            phi1 = -np.arccos(cos_phi1)

        return self.__class__(self.mu, r1, v1, phi1)
            
    def v_at_radius(self, r):
        return np.sqrt(2.0 * (self.energy + self.mu / r))

    @property
    def cos_nu(self):
        """Cosine of true anomaly"""
        return np.clip((self.p - self.r) / (self.e * self.r), -1.0, 1.0)

    @property
    def cos_E(self):
        """Cosine of eccentric anomaly"""
        cnu = self.cos_nu
        e = self.e

        return np.clip( (e + cnu) / (1.0 + e * cnu), -1.0, 1.0)
        
