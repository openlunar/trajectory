# using ForwardDiff, Plots, SparseArrays, MAT
# using LinearAlgebra, TrajectoryOptimization,
using StaticArrays
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
Base.size(::R3BP) = 12,3

struct CR3BP{T} <: AbstractModel
    μ::T # gravitational parameter for entire system
end
Base.size(::CR3BP) = 12,3

abstract type RK8 <: TrajectoryOptimization.Implicit
end

function discrete_dynamics(::Type{RK8}, model::AbstractModel, x::SVector{N,T}, u::SVector{M,T},t, dt::T) where {N,M,T}


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
    rm = x[ @SVector [7,8,9] ]
    vm = x[ @SVector [10,11,12] ]
    rdot = v
    vdot = -μ*r/norm(r)^3 - μ2*(r-rm)/norm(r-rm)^3 + u
    rmdot = vm
    vmdot = -(μ)*rm/norm(rm)^3
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
TrajectoryOptimization.∇²differential(::CR3BP, x, Qx) = Diagonal(@SVector zeros(6))


## Maneuver and discontinuity points
idx0 = 70000
idxf = 70742
findTraj(x,r_moon,v_moon,70000,idxf)

function findTraj(x,r_moon,v_moon,idx0,idxf)
    N = idxf - idx0
    n = 12 # number of states
    m = 3 # number of controls

    mu = 398600 #Gravitational Parameter
    mu2 = 4904.8695 #Moon's Gravitational Parameter
    model  = R3BP{Float64}(mu,mu2)

    U0 = [0.01*rand(m) for k = 1:N-1]; # initial control trajectory
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
    xtraj = [xx[k][1] for k = 1:N]
    ytraj = [xx[k][2] for k = 1:N]
    utraj = [norm(uu[k][:]) for k = 1:N-1]
    plotlyjs()
    plot(utraj)
    plot(xtraj,ytraj)
end

ΔV = 0.
for i = 3200:3499
    global utraj
    global ΔV += (utraj[i]+utraj[i-1])/2
end

###### Try using cubicSpline to downsample
xspline, tspline = cubicSpline(x,t,3001)

plot(xspline[1,:],xspline[2,:])

idx0 = 300
idxf = 1390
N = idxf - idx0
n = 6 # number of states
m = 3 # number of controls

mu = 398600 #Gravitational Parameter
model  = R2BP{Float64}(mu)

U0 = [0.01*rand(m) for k = 1:N-1]; # initial control trajectory
Q = 1.0*Diagonal(I,n)
Qf = 1.0*Diagonal(I,n)
R = 1.0e-1*Diagonal(I,m)

costfuns = map(1:N) do k #Brian helped me with this
    xf = xspline[:,idx0+k]
    LQRCost(Q,R,xf)
end

obj = Objective(costfuns)

#I don't currently have constraints 2/26/20
# constraints = Constraints(N) # define constraints at each time step
# for k = 1:N-1
    # constraints[k] += bnd
# end
# constraints[N] += goal

x0 = xspline[:,idx0+1]
xf = xspline[:,idxf]
t0 = tspline[idx0+1]
tf = tspline[idxf]
prob = Problem(model, obj, xf, tf, x0=x0, t0=t0) # construct problem
initial_controls!(prob,U0) # initialize problem with controls

# solver = solve!(prob, ALTROSolverOptions{Float64}())
solver = iLQRSolver(prob)

# solver = ALTROSolver(prob)
# solver = AugmentedLagrangianSolver(prob)

solve!(solver)

xx = states(solver)
uu = controls(solver)

N = prob.N
xtraj = [xx[k][1] for k = 1:N]
ytraj = [xx[k][2] for k = 1:N]
utraj = [norm(uu[k][:]) for k = 1:N-1]
plot(utraj)
plot(xtraj,ytraj)
