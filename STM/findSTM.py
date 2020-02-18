# from julia.api import Julia
# jl = Julia(compiled_modules=False)

from julia.api import Julia
api = LibJulia.load()
api.sysimage = "PATH/TO/CUSTOM/sys.so"
api.init_julia()

from julia import Main

from diffeqpy import de
import matplotlib.pyplot as plt
import numpy as np

def CR3BP(X, mu, t):
    Xdot = np.zeros(6)
    x,y,z,vx,vy,vz = X[:]

    mu1 = 1-mu
    r13 = ((x+mu)**2 + y**2 + z**2)**1.5 #distance from M1 to P cubed
    r23 = ((x-mu1)**2 + y**2 + z**2)**1.5 #distance from M2 to P cubed

    Xdot[0:3] = vx,vy,vz
    Xdot[3] = 2*vy + x - mu1*(x+mu)/r13 - mu*(x-mu1)/r23
    Xdot[4] = -2*vx + y - mu1*y/r13 - mu*y/r23
    Xdot[5] = -mu1*z/r13 - mu*z/r23
    if gradient:
        G = np.zeros((3,3))
        r15 = ((x+mu)**2 + y**2 + z**2)**2.5 #distance from M1 to P raised to the fifth power
        r25 = ((x-mu1)**2 + y**2 + z**2)**2.5 #distance from M2 to P raised to the fifth power

        G[0,0]      = 1 - mu1*(1/r13 - 3*(x+mu)**2/r15) - mu*(1/r23 - 3*(x-mu1)**2/r25);
        G[0,1]      = 3*mu1*(x+mu)*y/r15 + 3*mu*(x-mu1)*y/r25;
        G[0,2]      = 3*mu1*(x+mu)*z/r15 + 3*mu*(x-mu1)*z/r25;
        G[1,1]      = 1 - mu1*(1/r13 - 3*y**2/r15) - mu*(1/r23 - 3*y**2/r25);
        G[1,2]      = 3*mu1*y*z/r15 + 3*mu*y*z/r25;
        G[2,2]      = - mu1*(1/r13 - 3*z**2/r15) - mu*(1/r23 - 3*z**2/r25);

        G[1,0]      = G[0,1]; #Symmetric Matrix
        G[2,0]      = G[0,2];
        G[2,1]      = G[1,2];
        return Xdot, G
    else:
        return Xdot, None


def f(u,p,t):
    return -u

u0 = np.array([1,2])
tspan = (0., 10.)
prob = de.ODEProblem(f, u0, tspan)
sol = de.solve(prob)


plt.plot(sol.t,sol.u)
plt.show()
