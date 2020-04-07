dir = "/home/rexlab/trajectory/"
include(string(dir,"utils/OrbitDynamics.jl"))
include(string(dir,"utils/OptimizationSetup.jl"))

file = matopen(string(dir,"beresheet/BeresheetMoonSunEphem.mat"))
r_moon = read(file,"pos_moon")
v_moon = read(file,"vel_moon")
r_sun = read(file,"pos_sun")
v_sun = read(file,"vel_sun")
x = read(file,"x")
et = read(file,"times")
close(file)
t = et - et[1]*ones(size(et))

## Entire Trajectory
xx0, uu0 = findTraj(x,r_moon,v_moon,1,5317)
tout01, xx01 = ode8(R3BPdynamics,0,10500,xx0[end],h=60)
tout02, xx02 = ode8(R3BPdynamics,0,10500,xx01[end],h=60)
tout03, xx03 = ode8(R3BPdynamics,0,10500,xx02[end],h=60)
tout04, xx04 = ode8(R3BPdynamics,0,10500,xx03[end],h=60)
tout05, xx05 = ode8(R3BPdynamics,0,10500,xx04[end],h=60)
tout06, xx06 = ode8(R3BPdynamics,0,10500,xx05[end],h=60)
tout07, xx07 = ode8(R3BPdynamics,0,10500,xx06[end],h=60)
tout08, xx08 = ode8(R3BPdynamics,0,10500,xx07[end],h=60)
tout09, xx09 = ode8(R3BPdynamics,0,10500,xx08[end],h=60)

# scatter!([xx09[end-7][1],0],[xx09[end-7][2],0])
# scatter!([xx1[1][1],0],[xx1[1][2],0])

xx1, uu1 = findTraj(x,r_moon,v_moon,6892,36612)

# A = [xx0;xx01[2:end];xx02[2:end];xx03[2:end];xx04[2:end];xx05[2:end];xx06[2:end];xx07[2:end];xx08[2:end];xx09[2:end-7];xx1]
# A_array = convert2array(A)
# xxa, uua = findTraj(A_array[1:6,:],r_moon,v_moon,1,36604)
# xxa_array = convert2array(xxa)
# uua_array = convert2array(uua)
# plot(xxa_array[1,:],xxa_array[2,:])
# ua = [norm(uua[i][1:6]) for i = 1:length(uua)]
# plot(ua)

tout11, xx11 = ode8(R3BPdynamics,0,10500,xx1[end],h=60)
tout12, xx12 = ode8(R3BPdynamics,0,10500,xx11[end],h=60)
tout13, xx13 = ode8(R3BPdynamics,0,10500,xx12[end],h=60)
tout14, xx14 = ode8(R3BPdynamics,0,10500,xx13[end],h=60)
tout15, xx15 = ode8(R3BPdynamics,0,10500,xx14[end],h=60)
tout16, xx16 = ode8(R3BPdynamics,0,10500,xx15[end],h=60)

# B_array = convert2array(B)
# xxb, uub = findTraj(B_array[1:6,:],r_moon,v_moon,1,37641)
# ub = [norm(uub[i][1:6]) for i = 1:length(uub)]
# plot(ub)

# scatter!([B[end-13][1],0],[B[end-13][2],0])
# scatter!([xx2[1][1],0],[xx2[1][2],0])
# plot(xx11_array[1,:],xx11_array[2,:])

xx2, uu2 = findTraj(x,r_moon,v_moon,36613,70741)

xxtot = [xx0;xx01[2:end];xx02[2:end];xx03[2:end];xx04[2:end];xx05[2:end];xx06[2:end];xx07[2:end];xx08[2:end];xx09[2:end-7];xx1;xx11[2:end];xx12[2:end];xx13[2:end];xx14[2:end];xx15[2:end];xx16[2:end-15];xx2]

xxtot_array = convert2array(xxtot)
xxtot, uutot = findTraj(xxtot_array[1:6,1025:end],r_moon,v_moon,1,70742)

xxtot_array = convert2array(xxtot)
uutot_array = convert2array(uutot)
uutot_array[:,5600:6200] = zeros(6,601)
uutot_array[:,36400:36860] = zeros(6,461)
utot = [norm(uutot_array[1:3,i]) for i = 1:length(uutot)]

plot(xxtot_array[1,:],xxtot_array[2,:])
utot = [norm(uutot[i][1:6]) for i = 1:length(uutot)]
plot(utot,label="u")

#write to mat file
file = matopen("fullBeresheetTraj.mat", "w")
write(file, "x", xxtot_array[1:6,:])
write(file, "u", uutot_array[1:3,:])
write(file, "et", et)
write(file, "t", t)
close(file)

## write to spk file
file = matopen(string(dir,"beresheet/fullBeresheetTraj.mat"))
xx0 = read(file,"x")
uu0 = read(file,"u")
et = read(file,"et")
t = read(file,"t")
close(file)

spkfile = spkopn("/home/rexlab/trajectory/beresheet/fullBeresheetTraj.bsp", "Beresheet Trajectory SPK13 File", 1000)
spkw13(spkfile, -5441, 399, "J2000", et[2], et[end], "1", 3, convert2vector(xx0), et[2:end])
spkcls(spkfile)
furnsh(string(dir,"beresheet/fullBeresheetTraj.bsp"))
boddef("beresheet",-5441)
r1 = spkezr("beresheet",et[2],"J2000","none","earth")[1]
xxr = [spkezr("beresheet",idx,"J2000","none","earth")[1] for idx = et[2]:60:et[end]]

###### Try using cubicSpline to downsample
# xspline, tspline = cubicSpline(x,t,3001)
#
# plot(xspline[1,:],xspline[2,:])
#
# idx0 = 300
# idxf = 1390
# N = idxf - idx0
# n = 6 # number of states
# m = 3 # number of controls
#
# mu = 398600 #Gravitational Parameter
# model  = R2BP{Float64}(mu)
#
# U0 = [0.01*rand(m) for k = 1:N-1]; # initial control trajectory
# Q = 1.0*Diagonal(I,n)
# Qf = 1.0*Diagonal(I,n)
# R = 1.0e-1*Diagonal(I,m)
#
# costfuns = map(1:N) do k #Brian helped me with this
#     xf = xspline[:,idx0+k]
#     LQRCost(Q,R,xf)
# end
#
# obj = Objective(costfuns)
#
# #I don't currently have constraints 2/26/20
# # constraints = Constraints(N) # define constraints at each time step
# # for k = 1:N-1
#     # constraints[k] += bnd
# # end
# # constraints[N] += goal
#
# x0 = xspline[:,idx0+1]
# xf = xspline[:,idxf]
# t0 = tspline[idx0+1]
# tf = tspline[idxf]
# prob = Problem(model, obj, xf, tf, x0=x0, t0=t0) # construct problem
# initial_controls!(prob,U0) # initialize problem with controls
#
# # solver = solve!(prob, ALTROSolverOptions{Float64}())
# solver = iLQRSolver(prob)
#
# # solver = ALTROSolver(prob)
# # solver = AugmentedLagrangianSolver(prob)
#
# solve!(solver)
#
# xx = states(solver)
# uu = controls(solver)
#
# N = prob.N
# xtraj = [xx[k][1] for k = 1:N]
# ytraj = [xx[k][2] for k = 1:N]
# utraj = [norm(uu[k][:]) for k = 1:N-1]
# plot(utraj)
# plot(xtraj,ytraj)
