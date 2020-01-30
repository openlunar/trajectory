import numpy as np
from scipy.linalg import norm

from scipy.integrate import RK45

from propagate.forces import j2_gravity, gravity

def dynamics(t, x):
    f = np.zeros_like(x)

    #f[0:6] = j2_gravity(x[0:6], mu = 3.986004354360959e14, j2 = 1.08262668e-3, r_m = 6356800.0)
    #f[6:12] = j2_gravity(x[6:12], mu = 3.986004354360959e14, j2 = 1.08262668e-3, r_m = 6356800.0)
    f[0:6] = gravity(t, x[0:6], mu = 3.986004354360959e14, label='earth-sc')
    f[6:12] = gravity(t, x[6:12], mu = 3.986004354360959e14, label='earth-moon')

    f_moon = gravity(t, x[0:6] - x[6:12], mu = 4.902800066163796e12, label='moon-sc')
    f[3:6] += f_moon[3:6]
    return f

class Dynamics(object):
    # Constants: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/jgre.20097
    EARTH_J2_KWARGS = {'j2': 1.082626335439e-3, 'r_m': 6356800.0}
    MOON_J2_KWARGS  = {'j2': 2.0330530-4,       'r_m': 1738000.0}

    def __init__(self,
                 fun_earth = gravity,
                 fun_moon  = gravity,
                 earth_j2_kwargs = None,
                 moon_j2_kwargs  = None,
                 mu_earth = 3.986004354360959e14,
                 mu_moon  = 4.902800066163796e12,
                 j2_moon  = None,
                 j2_earth = None,
                 rm_moon  = None,
                 rm_earth = None):
        # Allow defaults to be overridden
        if j2_moon is not None:  self.MOON_J2_KWARGS['j2'] = j2_moon
        if j2_earth is not None: self.EARTH_J2_KWARGS['j2'] = j2_earth
        if rm_moon is not None:  self.MOON_J2_KWARGS['r_m'] = rm_moon
        if rm_earth is not None: self.EARTH_J2_KWARGS['r_m'] = rm_earth

        # Allow user to specify gravity functions to use between each pair of bodies
        if fun_earth == gravity:
            self.kwargs_earth = {}
        else:
            self.kwargs_earth = self.EARTH_J2_KWARGS

        if fun_moon == gravity:
            self.kwargs_moon = {}
        else:
            self.kwargs_moon  = self.MOON_J2_KWARGS
            
        self.fun_moon = fun_moon
        self.fun_earth = fun_earth
        self.mu_earth = mu_earth
        self.mu_moon = mu_moon

    def __call__(self, t, x):
        f = np.zeros_like(x)
        
        f[0:6] = self.fun_earth(t, x[0:6], mu = self.mu_earth,
                                **(self.kwargs_earth))
        f[6:12] = self.fun_earth(t, x[6:12], mu = self.mu_earth,
                                 **(self.kwargs_earth))

        f_moon = self.fun_moon(t, x[0:6] - x[6:12], mu = self.mu_moon,
                               **(self.kwargs_moon))
        
        f[3:6] += f_moon[3:6]

        return f
        
def propagate_to(fun, t0, x0, t1, plot = True, **kwargs):
    integ = RK45(fun, t0, x0, t1, **kwargs)

    mrs = []
    ers = []
    ts = []
    
    while integ.t < t1:
        if integ.status == 'failed':
            break
        else:
            integ.step()

        ts.append(integ.t)
        mrs.append(norm(integ.y[0:3] - integ.y[6:9]))
        ers.append(norm(integ.y[0:3]))

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(np.array(ts) - t0, mrs)
        ax.plot(np.array(ts) - t0, ers)
        ax.axhline(66183267.4546323)
        plt.show()

    return integ.y

def propagate_to_lunar_radius(fun, t0, x0, r2, t1_max,
                              max_step = 600.0,
                              min_step = 1e-4,
                              plot     = True,
                              **kwargs):
    integ = RK45(fun, t0, x0, t1_max, max_step = max_step, **kwargs)

   
    ts = []
    rs = []
    
    time   = t0
    radius = norm(integ.y[0:3] - integ.y[6:9])
    while radius > r2 and integ.status != 'failed':
        ts.append(time)
        rs.append(radius)

        try:
            integ.step()
        except ValueError:
            break
            
        prev_time   = time
        prev_radius = radius
        time        = integ.t
        radius      = norm(integ.y[0:3] - integ.y[6:9])

        denom = time - prev_time
        if denom < min_step:
            break
        rdot = (radius - prev_radius) / denom
        dt = -(radius - r2) / rdot
        print("dt = {}, step_size = {}".format(dt, integ.step_size))
        if dt > 0 and dt < max_step:
            integ.max_step = dt * 0.6

    ts.append(time)
    rs.append(radius)

    interp = integ.dense_output()
    print("dt = {}, step_size = {}".format(dt, integ.step_size))

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(np.array(ts) - t0, rs)
        ax.axhline(r2)

        x_soi = interp(integ.t + dt)
        ax.scatter([integ.t + dt - t0], [ norm(x_soi[0:3] - x_soi[6:9]) ], c='r')
        
        plt.show()
    
    
    return integ.t - dt, interp(integ.t - dt)

    

    
