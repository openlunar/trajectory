import unittest

from .context import *

class Test_PatchedConic(unittest.TestCase):

    def test_init_patched_conic(self):
        dlam1 = 0.01
        dv0 = 0.1
        leo = Orbit.circular(PatchedConic.mu_earth, 6371400.0 + 185000.0)

        opt_x0 = np.array([49.9 * np.pi/180.0, leo.v + 3140.0])
        opt_x1 = opt_x0 + np.array([dlam1, dv0])
        x0 = init_patched_conic(opt_x0, 1937000.0, leo.r, 0.0)
        x1 = init_patched_conic(opt_x1, 1937000.0, leo.r, 0.0)
        
        self.assertEqual(x0.depart.r, x1.depart.r)
        self.assertEqual(x0.depart.v + dv0, x1.depart.v)
        self.assertEqual(x0.lam1 + dlam1, x1.lam1)
        self.assertEqual(x0.rf, x1.rf)
        self.assertEqual(x0.depart.phi, x1.depart.phi)

        depart = Orbit(leo.mu, leo.r, leo.v + 3140.0)
        arrive = depart.at(x0.arrive.r, sign='+')
        x2 = PatchedConic(depart, arrive, lam1 = x0.lam1, rf = x0.rf)
        self.assertEqual(x0.depart.r, x2.depart.r)
        self.assertEqual(x0.depart.v, x2.depart.v)
        self.assertEqual(x0.lam1, x2.lam1)
        self.assertEqual(x0.rf, x2.rf)
        self.assertEqual(x0.depart.phi, x2.depart.phi)

    def test_dv0_gradients(self):
        dv0 = 0.1
        leo = Orbit.circular(PatchedConic.mu_earth, 6371400.0 + 185000.0)
        opt_x0 = np.array([49.9 * np.pi/180.0, leo.v + 3140.0])
        opt_x1 = opt_x0 + np.array([0.0, 0.1])
        x0 = init_patched_conic(opt_x0, 1937000.0, leo.r, 0.0)
        x1 = init_patched_conic(opt_x1, 1937000.0, leo.r, 0.0)
        x0g = PatchedConicGradients(x0)

        daf = x0g.daf_dv0 * dv0
        daf_ = x1.af - x0.af

        defm = x0g.def_dv0 * dv0
        defm_ = x1.ef - x0.ef

        dv1 = x0g.dv1_dv0 * dv0
        dv1_ = x1.arrive.v - x0.arrive.v

        drpl = x0g.drpl_dv0 * dv0
        drpl_ = x1.rpl - x0.rpl

        dvpl = x0g.dvpl_dv0 * dv0
        dvpl_ = x1.vpl - x0.vpl

        dphi1 = x0g.dphi1_dv0 * dv0
        dphi1_ = x1.arrive.phi - x0.arrive.phi

        dphi2 = x0g.dphi2_dv0 * dv0
        dphi2_ = x1.phi - x0.phi

        dg = x0g.dg_dv0 * dv0
        dg_ = x1.g - x0.g

        df = x0g.df_dv0 * dv0
        df_ = x1.f - x0.f

        deltav1 = x0g.ddeltav1_dv0 * dv0
        deltav1_ = x1.deltav1 - x0.deltav1

        dv2 = x0g.dv2_dv0 * dv0
        dv2_ = x1.v - x0.v

        np.testing.assert_approx_equal(dv2, dv2_, significant=3)
        np.testing.assert_approx_equal(deltav1, deltav1_, significant=3)
        np.testing.assert_approx_equal(dphi2, dphi2_, significant=3)
        np.testing.assert_approx_equal(dphi1, dphi1_, significant=2)
        np.testing.assert_approx_equal(dvpl, dvpl_, significant=2)
        np.testing.assert_approx_equal(drpl, drpl_, significant=2)
        np.testing.assert_approx_equal(dv1, dv1_, significant=3)
        np.testing.assert_approx_equal(defm, defm_, significant=3)
        np.testing.assert_approx_equal(daf, daf_, significant=3)
        np.testing.assert_approx_equal(dg, dg_, significant=2)
        np.testing.assert_approx_equal(df, df_, significant=2)
