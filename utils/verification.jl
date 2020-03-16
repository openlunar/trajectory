# include("/home/rexlab/trajectory/utils/OrbitDynamics.jl")
include("OrbitDynamics.jl")

## verify rk8() and ode8()
# file = matopen("/home/rexlab/trajectory/beresheet/BeresheetMoonSunEphem.mat")
file = matopen("beresheet/BeresheetMoonSunEphem.mat")
r_moon = read(file,"pos_moon")
v_moon = read(file,"vel_moon")
r_sun = read(file,"pos_sun")
v_sun = read(file,"vel_sun")
x = read(file,"x")
t = read(file,"times")
close(file)

# rv, theta = J2000_2_CR3BP(x,t,[r_moon;v_moon])

# rv0 = rv[:,1]
# t0 = theta[1]
# tf = 5
# rk8(CR3BPdynamics,rv[:,1],theta[1],theta[2]-theta[1])
# tout, yout = ode8(CR3BPdynamics,t0,tf,rv0)

# plot(yout[:,1],yout[:,2])
# plot!(rv[1,:],rv[2,:])

idx = [1 9672+200 19349+200 36626+200 60000+200]
t = t - t[1]*ones(size(t))

tt1, xx1 = ode8(R3BPdynamics,t[idx[1]],t[idx[2]],[x[:,idx[1]];r_moon[:,idx[1]];v_moon[:,idx[1]]])
tt2, xx2 = ode8(R3BPdynamics,t[idx[2]],t[idx[3]],[x[:,idx[2]];r_moon[:,idx[2]];v_moon[:,idx[2]]])
tt3, xx3 = ode8(R3BPdynamics,t[idx[3]],t[idx[4]],[x[:,idx[3]];r_moon[:,idx[3]];v_moon[:,idx[3]]])
tt4, xx4 = ode8(R3BPdynamics,t[idx[4]],t[idx[5]],[x[:,idx[4]];r_moon[:,idx[4]];v_moon[:,idx[4]]])
tt5, xx5 = ode8(R3BPdynamics,t[idx[5]],t[end],   [x[:,idx[5]];r_moon[:,idx[5]];v_moon[:,idx[5]]])

tt, xx = ode8(R3BPdynamics,t[36626],t[19349],[x[:,36626];r_moon[:,36626];v_moon[:,36626]])

plot(x[1,:],x[2,:])
plot(xx1[:,1],xx1[:,2])
plot!(xx2[:,1],xx2[:,2])
plot!(xx3[:,1],xx3[:,2])
plot!(xx4[:,1],xx4[:,2])
plot!(xx5[:,1],xx5[:,2])
