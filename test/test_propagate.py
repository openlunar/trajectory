import unittest

from .context import *

class Test_Propagate(unittest.TestCase):

    def setUp(self):
        
        loader = SpiceLoader()

        # Arbitrary arrival date at SOI (assuming model)
        arrival_date = '21JUN2022'
        arrival_time = spice.utc2et(arrival_date)
        
        init = InitialState(arrival_time)
        print("Initial ECI state is {}".format(init.x_depart_pre))


        dynamics = Dynamics(fun_earth = j2_gravity,
                            fun_moon  = gravity)
        x0 = np.hstack((init.x_depart_post, init.x_moon_depart, np.identity(6).reshape(36)))
    
        t, x, Phi = propagate_to_lunar_radius(dynamics, init.depart_time, np.array(x0),
                                              PatchedConic.r_soi, arrival_time + 2 * 24 * 3600.0, plot = True)

        self.dynamics = dynamics
        self.x0 = x0
        self.t0 = init.depart_time
        self.t = t
        self.x = x
        self.Phi = Phi

    def test_propagate_to(self):
        """Time-based propagation should give the same final state as to-lunar-radius propagation"""
        t, x, Phi = propagate_to(self.dynamics, self.t0, np.array(self.x0), self.t, plot = False)
        np.testing.assert_allclose(x[0:6], self.x[0:6], atol=1e-1)

    def test_propagate_to_lunar_radius(self):
        """Lunar-radius-based propagation should take us to the fixed lunar sphere of influence"""
        np.testing.assert_approx_equal(norm(self.x[0:3] - self.x[6:9]), PatchedConic.r_soi, significant=4)

    def test_Phi(self):
        
        for ii in range(0,6):
            x0 = np.array(self.x0)
            x0[ii] += 1.0

            # How much x0 deviates from the original:
            dx0 = x0[0:6] - self.x0[0:6]
            
            t, x, Phi = propagate_to(self.dynamics, self.t0, x0, self.t, plot = True)

            stm_dx = Phi.dot(dx0)
            true_dx = x[0:6] - self.x[0:6]

            import matplotlib.pyplot as plt
            plt.show()

            import pdb
            pdb.set_trace()

            np.testing.assert_allclose(stm_dx, true_dx)
            
