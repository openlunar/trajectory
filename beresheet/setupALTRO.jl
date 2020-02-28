using ForwardDiff, Plots, SparseArrays, MAT
using LinearAlgebra, TrajectoryOptimization, StaticArrays
include("cubicSpline.jl")

function readmat(filename)
    cd("/home/rexlab/trajectory/beresheet/")
    file = matopen(filename)
    xt = read(file,"BeresheetTraj")
    close(file)
    N,M = size(xt)
    x = reshape(xt[1:N,1:6]',6,N)
    if M == 6
        return x
    elseif M == 7
        t = xt[1:N,7]*86400 #sec
        return x,t
    else
        error(".mat file data should be Nx6 or Nx7")
    end
end

struct R2BP{T} <: AbstractModel
    mu::T #gravitational parameter
end
Base.size(::R2BP) = 6,3

function TrajectoryOptimization.dynamics(model::R2BP, x, u) # inplace dynamics
    μ = model.mu
    r = x[ @SVector [1,2,3] ]
    v = x[ @SVector [4,5,6] ]
    rdot = v
    vdot = -μ*r/norm(r)^3 + u
    return [rdot; vdot]
end

TrajectoryOptimization.∇²differential(::R2BP, x, Qx) = Diagonal(@SVector zeros(6))

x, t = readmat("BeresheetTrajectory.mat")
idx0 = 5000
idxf = 8000
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
    xf = x[:,idx0+k]
    LQRCost(Q,R,xf)
end

obj = Objective(costfuns)

#I don't currently have constraints 2/26/20
# constraints = Constraints(N) # define constraints at each time step
# for k = 1:N-1
    # constraints[k] += bnd
# end
# constraints[N] += goal

x0 = x[:,idx0+1]
xf = x[:,idxf]
t0 = t[idx0+1]
tf = t[idxf]
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
