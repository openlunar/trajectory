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


        dynamics = Dynamics(fun_earth = gravity,
                            fun_moon  = gravity)
        x0 = np.hstack((init.x_depart_post, init.x_moon_depart, np.identity(6).reshape(36)))
    
        t, x, Phi = propagate_to_lunar_radius(dynamics, init.depart_time, np.array(x0),
                                              PatchedConic.r_soi, arrival_time + 2 * 24 * 3600.0, plot = True,
                                              method = 'DOP853')

        self.dynamics = dynamics
        self.x0 = x0
        self.t0 = init.depart_time
        self.t = t
        self.x = x
        self.Phi = Phi

        self.rtol = 1e-9
        self.atol = np.ones(48) * 1e-15
        self.atol[3:6] = 1e-18
        self.atol[12:] = 1e-18

    def test_propagate_to(self):
        """Time-based propagation should give the same final state as to-lunar-radius propagation"""
        t, x, Phi = propagate_to(self.dynamics, self.t0, np.array(self.x0), self.t, plot = False,
                                 integrator = spint.DOP853, rtol=self.rtol, atol=self.atol)
        np.testing.assert_allclose(x[0:6], self.x[0:6], atol=1e-1)

    def test_propagate_to_lunar_radius(self):
        """Lunar-radius-based propagation should take us to the fixed lunar sphere of influence"""
        np.testing.assert_approx_equal(norm(self.x[0:3] - self.x[6:9]), PatchedConic.r_soi, significant=4)

    def test_Phi_single_step(self):
        """Propagated state transition model should agree with single-step state propagation"""
        Phi_num = np.zeros((6,6))

        for ii in range(-1,6):
            x0 = np.array(self.x0)

            if ii >= 0:
                x0[ii] += 1.0

            # How much x0 deviates from the original:
            dx0 = x0[0:6] - self.x0[0:6]
            print("dx0 = {}".format(dx0))

            integ = spint.DOP853(self.dynamics, self.t0, x0, self.t0 + 1.0, rtol=self.rtol, atol=self.atol)
            integ.step()

            if ii < 0:
                ref_x1 = np.array(integ.y[0:6])
                Phi = integ.y[12:].reshape((6,6))
            else:
                print("dx1 = {}".format(integ.y[0:6] - ref_x1))
                Phi_num[:,ii] = integ.y[0:6] - ref_x1

        np.testing.assert_allclose(Phi, Phi_num, atol=1e-12, rtol=0.5)
        
    def test_Phi(self):
        """Propagated state transition model should agree with state propagation"""
        Phi_num = np.zeros((6,6))
        
        for ii in range(0,6):
            x0 = np.array(self.x0)
            x0[ii] += 0.1

            # How much x0 deviates from the original:
            dx0 = x0[0:6] - self.x0[0:6]
            
            t, x, Phi = propagate_to(self.dynamics, self.t0, x0, self.t,
                                     label = str(ii), rtol=self.rtol, atol=self.atol)

            stm_dx = self.Phi.dot(dx0)
            true_dx = x[0:6] - self.x[0:6]

            Phi_num[:,ii] = true_dx / 0.1

            np.testing.assert_allclose(stm_dx, true_dx, rtol=5e-4, atol=0.1)

