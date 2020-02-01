import numpy as np
from scipy.linalg import norm

from scipy.integrate import RK45

from propagate.forces import j2_gravity, gravity


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
            G += G_moon
            
            # Compute state transition matrix time derivative
            Phi = x[12:48].reshape((6,6))
            F   = np.vstack(( np.hstack(( np.zeros((3,3)), np.identity(3)  )),
                              np.hstack(( G,               np.zeros((3,3)) )) ))

            Phi_dot  = F.dot(Phi)
            f[12:48] = Phi_dot.reshape(36)

        return f
        
def propagate_to(fun, t0, x0, t1, plot = True, **kwargs):
    integ = RK45(fun, t0, x0, t1, **kwargs)

    mrs = []
    ers = []
    ts = []
    
    while integ.t < t1:
        if integ.status == 'failed':
            raise ValueError("integration failed")
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

    return integ.t, integ.y[0:12], integ.y[12:].reshape((6,6))

def propagate_to_lunar_radius(fun, t0, x0, r2, t1_max,
                              max_step = 600.0,
                              min_step = 1e-5,
                              plot     = True,
                              axes     = None,
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
            print("Warning: Integration step failure")
            break
            
        prev_time   = time
        prev_radius = radius
        time        = integ.t
        radius      = norm(integ.y[0:3] - integ.y[6:9])

        num   = radius - prev_radius
        denom = time - prev_time
        if denom < min_step:
            break
        rdot = num / denom
        dt = -(radius - r2) / rdot
        
        print("dt = {}, step_size = {}".format(dt, integ.step_size))
        if dt > 0 and dt < max_step:
            integ.max_step = dt

    ts.append(time)
    rs.append(radius)

    interp = integ.dense_output()
    print("dt = {}, step_size = {}".format(dt, integ.step_size))

    if plot:
        import matplotlib.pyplot as plt
        if not axes:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        axes.plot(np.array(ts) - t0, rs)
        axes.axhline(r2)

        x_soi = interp(integ.t - dt)
        axes.scatter([integ.t - dt - t0], [ norm(x_soi[0:3] - x_soi[6:9]) ], c='r')

        x_2 = interp(integ.t + dt)
        axes.scatter([integ.t + dt - t0], [ norm(x_2[0:3] - x_2[6:9]) ], c='b')
        
        #plt.show()
    

    y = interp(integ.t - dt)
    return integ.t - dt, y[0:12], y[12:].reshape((6,6))

    

    
