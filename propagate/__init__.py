import numpy as np
from scipy.linalg import norm

from scipy.integrate import RK45, DOP853, solve_ivp

from propagate.forces import j2_gravity, gravity, zero_gravity


class Dynamics(object):
    # Constants: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/jgre.20097

    def __init__(self,
                 fun_earth = gravity,
                 fun_moon  = gravity,
                 mu_earth  = 3.986004354360959e14,
                 mu_moon   = 4.902800066163796e12,
                 j2_moon   = 2.0330530-4,
                 j2_earth  = 1.082626335439e-3,
                 req_moon  = 1738000.0,
                 req_earth = 6356800.0):
        # Allow defaults to be overridden
        self.kwargs_earth = {'j2': j2_earth, 'r_eq': req_earth }
        self.kwargs_moon =  {'j2': j2_moon,  'r_eq': req_moon }
            
        self.fun_moon = fun_moon
        self.fun_earth = fun_earth
        self.mu_earth = mu_earth
        self.mu_moon = mu_moon

    def __call__(self, t, x):
        f = np.zeros_like(x)

        gradient = True if len(f) > 12 else False
        
        f[0:6], G = self.fun_earth(x[0:6], self.mu_earth,
                                   gradient = gradient,
                                   **(self.kwargs_earth))
        f[6:12], tmp = self.fun_earth(x[6:12], self.mu_earth,
                                 gradient = False, # No STM needed for lunar motion
                                 **(self.kwargs_earth))

        f_moon, G_moon = self.fun_moon(x[0:6] - x[6:12], self.mu_moon,
                                       gradient = gradient,
                                       **(self.kwargs_moon))
        
        f[3:6] += f_moon[3:6]

        if len(f) > 12:
            if self.fun_moon != zero_gravity:
                G += G_moon
            
            # Compute state transition matrix time derivative
            Phi = x[12:48].reshape((6,6))
            F   = np.vstack(( np.hstack(( np.zeros((3,3)), np.identity(3)  )),
                              np.hstack(( G,               np.zeros((3,3)) )) ))

            Phi_dot  = F.dot(Phi)
            f[12:48] = Phi_dot.reshape(36)

        return f

def default_atol():
    atol = np.ones(48) * 1e-15
    atol[3:6] *= 1e-3
    return atol
    
def propagate_to(fun, t0, x0, t1, plot = False, integrator = DOP853, axes = None, label = '',
                 atol = None,
                 **kwargs):
    if atol is None:
        atol = default_atol()
    
    integ = integrator(fun, t0, x0, t1, **kwargs)
    
    mrs = []
    ers = []
    ts = []

    count = 0
    while integ.t < t1:
        if integ.status == 'failed':
            raise ValueError("integration failed")
        else:
            integ.step()
            count += 1

        ts.append(integ.t)
        mrs.append(norm(integ.y[0:3] - integ.y[6:9]))
        ers.append(integ.y[0:3])

    print("integration step count: {}".format(count))
    if plot:
        ers = np.vstack(ers).T
        
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111, projection='3d')
        #axes.plot(np.array(ts) - t0, mrs)
        axes.plot(ers[0,:], ers[1,:], ers[2,:], label = label, alpha=0.6)
        axes.axhline(66183267.4546323)

    return integ.t, integ.y[0:12], integ.y[12:].reshape((6,6))

def plot_solve_ivp_result(result, axes = None):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    if not axes:
        fig = plt.figure()
        axes = fig.add_subplot(111, projection='3d')

    y = result.y[0:3,:]# - result.y[6:9,:]
    m = result.y[6:9,:]
    
    axes.plot(y[0,:], y[1,:], y[2,:], label='trajectory')
    axes.plot(m[0,:], m[1,:], m[2,:], label='moon')

    mx, my, mz = result.y[6:9,-1]
    print("Moon pos: {}, {}, {}".format(mx, my, mz))
    axes.plot([0, mx], [0, my], [0, mz], label='arrival')
    return axes

def propagate_to_lunar_radius(fun, t0, x0, r2, t1_max,
                              max_step   = 600.0,
                              plot       = False,
                              axes       = None,
                              method     = 'DOP853',
                              rtol       = 1e-9,
                              atol       = None):

    def lunar_radius(t, y):
        return norm(y[0:3] - y[6:9]) - r2
    lunar_radius.terminal = True

    if atol is None:
        atol = default_atol()
    
    result = solve_ivp(fun, (t0, t1_max), x0, method = method, events = lunar_radius,
                       rtol=rtol, atol=atol, max_step=max_step)

    if plot:
        axes = plot_solve_ivp_result(result, axes)
        axes.axhline(r2)       

    return result.t[-1], result.y[0:12,-1], result.y[12:,-1].reshape((6,6))

def propagate_to_periselene(fun, t0, x0,
                            t1_max   = None,
                            max_step = 600.0,
                            method   = 'DOP853',
                            plot     = False,
                            axes     = None,
                            rtol     = 1e-9,
                            atol     = None):

    def periselene(t, y):
        x_lci = y[0:6] - y[6:12]
        r_hat = x_lci[0:3] / norm(x_lci[0:3])
        v_hat = x_lci[3:6] / norm(x_lci[3:6])
        return np.dot(r_hat, v_hat) < 1e-6
    periselene.terminal = True
    
    if atol is None:
        atol = default_atol()
    
    result = solve_ivp(fun, (t0, t1_max), x0, method = method, events = periselene,
                       rtol=rtol, atol=atol, max_step=max_step)

    if plot:
        axes = plot_solve_ivp_result(result, axes)
    
    return result.t[-1], result.y[0:12,-1], result.y[12:,-1].reshape((6,6))
