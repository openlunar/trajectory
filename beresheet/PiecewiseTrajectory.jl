# using ForwardDiff, Plots, SparseArrays, MAT
# using LinearAlgebra, TrajectoryOptimization,
using StaticArrays
using DataFrames, CSV
using TrajectoryOptimization
import TrajectoryOptimization.discrete_dynamics
include("/home/rexlab/trajectory/utils/OrbitDynamics.jl")
# include("setupALTRO.jl")

file = matopen("/home/rexlab/trajectory/beresheet/BeresheetMoonSunEphem.mat")
r_moon = read(file,"pos_moon")
v_moon = read(file,"vel_moon")
r_sun = read(file,"pos_sun")
v_sun = read(file,"vel_sun")
x = read(file,"x")
t = read(file,"times")
close(file)

t = t - t[1]*ones(size(t))

struct R3BP{T} <: AbstractModel
    μ::T #Earth's gravitational parameter
    μ2::T #Moon's gravitational parameter
end
Base.size(::R3BP) = 12,6

struct R4BP{T} <: AbstractModel
    μ::T #Earth's gravitational parameter
    μ2::T #Moon's gravitational parameter
    μ3::T #Sun's gravitational parameter
end
Base.size(::R4BP) = 18,9

struct CR3BP{T} <: AbstractModel
    μ::T # gravitational parameter for entire system
end
Base.size(::CR3BP) = 6,3

abstract type RK8 <: TrajectoryOptimization.Implicit
end


function discrete_dynamics(::Type{RK3}, model::AbstractModel, y::SVector{N,T}, u::SVector{M,T},t, dt::T) where {N,M,T}
    α    = @SVector [ 2/27, 1/9, 1/6, 5/12, .5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ];
    β    = @SMatrix [  2/27       0       0      0        0         0       0         0     0      0     0 0 0;
    1/36       1/12    0      0        0         0       0         0     0      0     0 0 0;
    1/24       0       1/8    0        0         0       0         0     0      0     0 0 0;
    5/12       0       -25/16 25/16    0         0       0         0     0      0     0 0 0;
    .05        0       0      .25      .2        0       0         0     0      0     0 0 0;
    -25/108    0       0      125/108  -65/27    125/54  0         0     0      0     0 0 0;
    31/300     0       0      0        61/225    -2/9    13/900    0     0      0     0 0 0;
    2          0       0      -53/6    704/45    -107/9  67/90     3     0      0     0 0 0;
    -91/108    0       0      23/108   -976/135  311/54  -19/60    17/6  -1/12  0     0 0 0;
    2383/4100  0       0      -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0 0;
    3/205      0       0      0        0         -6/41   -3/205    -3/41 3/41   6/41  0 0 0;
    -1777/4100 0       0      -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1 0]
    χ     = @SVector[0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840 ]; # Difference between two bottom layers of butcher tableau

    f = y*zeros(1,13)

    f[:,1] = dynamics(model,y,u,t)
    for j = 1:12
        f[:,j+1] = dynamics(model, y + dt*f*β[j,:],u,t + α[j]*dt)
        # f[:,j+1] = f[:,1]  - f[:,1]
    end
    return y + dt*f*χ
    # y-y

    # k1 = dynamics(model, x,             u, t       )*dt;
    # k2 = dynamics(model, x + k1/2,      u, t + dt/2)*dt;
    # k3 = dynamics(model, x - k1 + 2*k2, u, t       )*dt;
    # x + (k1 + 4*k2 + k3)/6
    # myvector

end

function TrajectoryOptimization.dynamics(model::R3BP, x, u, t) # inplace dynamics
    μ = model.μ
    μ2 = model.μ2
    r = x[ @SVector [1,2,3] ]
    v = x[ @SVector [4,5,6] ]
    u1 = u[ @SVector [1,2,3] ]
    u2 = u[ @SVector [4,5,6] ]
    rm = x[ @SVector [7,8,9] ]
    vm = x[ @SVector [10,11,12] ]
    rdot = v
    vdot = -μ*r/norm(r)^3 - μ2*(r-rm)/norm(r-rm)^3 + u1
    rmdot = vm
    vmdot = -(μ)*rm/norm(rm)^3 + u2
    return [rdot; vdot; rmdot; vmdot]
end

# function TrajectoryOptimization.dynamics(model::R3BP, x, u) # inplace dynamics
#     μ = model.μ
#     μ2 = model.μ2
#     r = x[ @SVector [1,2,3] ]
#     v = x[ @SVector [4,5,6] ]
#     rm = x[ @SVector [7,8,9] ]
#     vm = x[ @SVector [10,11,12] ]
#     rdot = v
#     vdot = -μ*r/norm(r)^3 - μ2*(r-rm)/norm(r-rm)^3 + u
#     rmdot = vm
#     vmdot = -(μ)*rm/norm(rm)^3
#     return [rdot; vdot; rmdot; vmdot]
# end

# function TrajectoryOptimization.dynamics(model::CR3BP, x, u) # inplace dynamics
#     μ = model.μf[:,j+1] = f[:,1]
#     return CR3BPdynamics(0,x,mu)
# end

TrajectoryOptimization.∇²differential(::R3BP, x, Qx) = Diagonal(@SVector zeros(12))
# TrajectoryOptimization.∇²differential(::CR3BP, x, Qx) = Diagonal(@SVector zeros(6))


## Maneuver and discontinuity points
# idx0 = 70000
# idxf = 70742
# findTraj(x,r_moon,v_moon,idx0,idxf)

# idx = [1, 9672, 19349, 36626]


function findTraj(x,r_moon,v_moon,idx0,idxf, U0 = [0.01*rand(6) for k = 1:idxf-idx0-1])
    N = idxf - idx0
    n = 12 # number of states
    m = 6 # number of controls

    model  = R3BP{Float64}(earth.μ, moon.μ)
    # model  = R4BP{Float64}(earth.μ, moon.μ, sun.μ)

    # U0 = [0.01*rand(m) for k = 1:N-1]; # initial control trajectory
    Q = 1.0*Diagonal(I,n)
    Qf = 1.0*Diagonal(I,n)
    R = 1.0e-1*Diagonal(I,m)

    costfuns = map(1:N) do k #Brian helped me with this
        xf = x[:,idx0+k]
        rfm = r_moon[:,idx0+k]
        vfm = v_moon[:,idx0+k]
        LQRCost(Q,R,[xf;rfm;vfm])
    end

    obj = Objective(costfuns)

    #I don't currently have constraints 2/26/20
    # constraints = Constraints(N) # define constraints at each time step
    # for k = 1:N-1
        # constraints[k] += bnd
    # end
    # constraints[N] += goal

    x0 = [x[:,idx0+1]; r_moon[:,idx0+1]; v_moon[:,idx0+1]]
    xf = [x[:,idxf]; r_moon[:,idxf]; v_moon[:,idxf]]
    t0 = t[idx0+1]
    tf = t[idxf]
    prob = Problem(model, obj, xf, tf, x0=x0, t0=t0) # construct problem
    initial_controls!(prob,U0) # initialize problem with controls

    # solver = solve!(prob, ALTROSolverOptions{Float64}())
    solver = iLQRSolver(prob)

    # solver = ALTROSolver(prob)
    # solver = AugmentedLagrangianSolver(prob)

    TrajectoryOptimization.solve!(solver)

    xx = states(solver)
    uu = controls(solver)

    N = prob.N
    # xtraj = [xx[k][1] for k = 1:N]
    # ytraj = [xx[k][2] for k = 1:N]
    # ztraj = [xx[k][3] for k = 1:N]
    # xmtraj = [xx[k][7] for k = 1:N]
    # ymtraj = [xx[k][8] for k = 1:N]
    # zmtraj = [xx[k][9] for k = 1:N]
    # u1traj = [norm(uu[k][1:3]) for k = 1:N-1]
    # u2traj = [norm(uu[k][4:6]) for k = 1:N-1]
    # return xtraj, ytraj, ztraj, xmtraj, ymtraj, zmtraj, u1traj, u2traj,
    return xx, uu
    # plotlyjs()
    # plot(utraj)
    # plot!(xtraj,ytraj)
end

function convert2array(xx)
    N = length(xx)
    M = length(xx[1])
    xx_array = zeros(M,N)
    for i = 1:M
        xx_array[i,:] = [xx[k][i] for k = 1:N]
    end
    return xx_array
end

function convert2vector(xx_array)
    m,n = size(xx_array)
    xx_vector = [xx_array[:,i] for i = 1:n]
end

# xx, uu = findTraj(x,r_moon,v_moon,59000,61000)
# xx_array = convert2array(xx)
# plot(xx_array[1,:],xx_array[2,:],xx_array[3,:])
#
# u1traj = [norm(uu[i][1:3]) for i = 1:length(uu)]
# u2traj = [norm(uu[i][4:6]) for i = 1:length(uu)]
#
#
# tt = (1:1999)/60
# plot(tt,u1traj,label='u',xlabel="time from epoch [hr]", ylabel="control magnitude km/s^2")
# plot!(tt,u2traj,label='u',xlabel="time from epoch [hr]", ylabel="control magnitude km/s^2")
#
# Δv = 0
# UU = 0
# for n = 1:1500
#     global UU += (uu[n+1] + uu[n])*60/2
# end
# Δv = UU[1:3]

## Entire Trajectory

function plotregion(idx0,idxf)
    xx, uu = findTraj(x,r_moon,v_moon,idx0,idxf)
    xx_array = convert2array(xx)
    u = [norm(uu[i][1:3]) for i = 1:length(uu)]
    uu_array = convert2array(uu)
    tt = (1:length(xx))/60
    p1 = plot(xx_array[1,:],xx_array[2,:],label="trajectory")
    # p2 = plot(u,label='u',xlabel="time from epoch [hr]", ylabel="control magnitude km/s^2")
    # plot(p1, p2, layout = (2, 1), legend = false)
end

x_norm = [norm(x[1:3,i]) for i = 1:70742]
# plot(x_norm)

function plot_chain(x_chain)
    n = length(x_chain)
    xtot = x_chain[1]
    if n > 1
        for i = 2:n
            xtot = [xtot; x_chain[i][2:end]]
        end
    end
    xx_array = convert2array(xtot)
    plot(xx_array[1,:],xx_array[2,:])
end

function make_chain(x_chain)
    n = length(x_chain)
    xtot = x_chain[1]
    if n > 1
        for i = 2:n
            xtot = [xtot; x_chain[i][2:end]]
        end
    end
    return xtot
end

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

write("xxtot_array.out",xxtot_array)
write("xxtot.out",xxtot)
write("uutot_array.out",uutot_array)
write("utot.out",utot)

CSV.write("xxtot_array.csv",  DataFrame(xxtot_array), writeheader=false)
CSV.write("uutot_array.csv",  DataFrame(uutot_array), writeheader=false)

#
# C = [xx16[1:end-15];xx2]
# C_array = [zeros(12,70742-34304+15) convert2array(C)]
# xxc, uuc = findTraj(C_array[1:6,:],r_moon,v_moon,70742-34304+16,70742)
#
# xxc_array = convert2array(xxc)
# uuc_array = convert2array(uuc)
# plot(xxc_array[1,:],xxc_array[2,:])
# uc = [norm(uuc[i][1:6]) for i = 1:length(uuc)]
# plot(uc)
#
# xx0_array = convert2array(xx0)
# uu0_array = convert2array(uu0)
# plot(xx0_array[1,:],xx0_array[2,:])
# u0 = [norm(uu0[i][1:6]) for i = 1:length(uu0)]
# plot(u0)
#
# xx1_array = convert2array(xx1)
# uu1_array = convert2array(uu1)
# plot(xx1_array[1,:],xx1_array[2,:])
# u1 = [norm(uu1[i][1:6]) for i = 1:length(uu1)]
# plot(u1)
#
# xx2_array = convert2array(xx2)
# uu2_array = convert2array(uu2)
# plot(xx2_array[1,:],xx2_array[2,:])
# u2 = [norm(uu2[i][1:6]) for i = 1:length(uu2)]
# plot(u2)
#
# p1 = plot(xx_array[1,:],xx_array[2,:]) # Make a line plot
# p2 = scatter(xx[1][1], xx[1][2]) # Make a scatter plot
# p3 = plot(x, y, xlabel = "This one is labelled", lw = 3, title = "Subtitle")
# p4 = histogram(x, y) # Four histograms each with 10 points? Why not!
# plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
#
# uu_tot = uu_array
# for i = 56:69
#     xx, uu = findTraj(x,r_moon,v_moon,1000*i,1000*(i+1))
#     xx_array = convert2array(xx)
#     global xx_tot = [xx_tot xx_array]
#     global uu_tot = [uu_tot uu_array]
# end
# xx, uu = findTraj(x,r_moon,v_moon,70000,70742)
# xx_array = convert2array(xx)
# xx_tot = [xx_tot xx_array]
# uu_tot = [uu_tot uu_array]

plot(xx_tot[1,:],xx_tot[2,:])
uu_tot_norm = [norm(uu_tot[4:6,k]) for k = 1:size(uu_tot)[2]]
plot(uu_tot_norm)

U0 = convert2vector(uu_tot)
xx_big, uu_big = findTraj(xx_tot[1:6,:],xx_tot[7:9,:],xx_tot[10:12,:],1,10000,U0)
xx_big = convert2array(xx_big)
uu_big = convert2array(uu_big)
plot(xx_big[1,:],xx_big[2,:])
plot!(xx_tot[1,:],xx_tot[2,:])
plot(uu_big[1,:],uu_big[2,:])

help

#
# ΔV = 0.
# for i = 3200:3499
#     global utraj
#     global ΔV += (utraj[i]+utraj[i-1])/2
# end

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
