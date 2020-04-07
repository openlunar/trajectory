struct EPHEM{T} <: AbstractModel
    μ::T #Earth's gravitational parameter
    μ2::T #Moon's gravitational parameter
    μ3::T #Sun's gravitational parameter
    r_moon::ScaledInterpolation{Array{Float64,1},1,Interpolations.BSplineInterpolation{Array{Float64,1},1,Array{Array{Float64,1},1},BSpline{Linear},Tuple{Base.OneTo{Int64}}},BSpline{Linear},Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}}
    R::ScaledInterpolation{Array{Float64,2},1,Interpolations.BSplineInterpolation{Array{Float64,2},1,Array{Array{Float64,2},1},BSpline{Linear},Tuple{Base.OneTo{Int64}}},BSpline{Linear},Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}}
end
Base.size(::EPHEM) = 6,3

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

#Runge-Kutta 8th order integration scheme
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
    end
    return y + dt*f*χ
end

# Propagating the position of the moon along with the position of the S.C.
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

function TrajectoryOptimization.dynamics(model::EPHEM, x, u, t) # inplace dynamics
    μ = model.μ     #Earth gravitational parameter
    μ2 = model.μ2   #Moon gravitational parameter

    # r = x[ @SVector [1,2,3] ]   #S.C. position
    # v = x[ @SVector [4,5,6] ]   #S.C. velocity
    # u = u[ @SVector [1,2,3] ]   #control input

    r = x[1:3]   #S.C. position
    v = x[4:6]   #S.C. velocity
    u = u[1:3]   #control input

    r2 = model.r_moon(t)
    R = model.R(t)

    rdot = v
    # vdot = -μ*r/norm(r)^3 - accel_gravity(1000*(r-r2),Array(R),20,20) + u #Using SatelliteDynamics.jl gravity model
    vdot = -μ*r/norm(r)^3 - μ2*(r-r2)/norm(r-r2)^3 + u #Using point mass model
    return [rdot; vdot]
end

TrajectoryOptimization.∇²differential(::R3BP, x, Qx) = Diagonal(@SVector zeros(12))
TrajectoryOptimization.∇²differential(::EPHEM, x, Qx) = Diagonal(@SVector zeros(6))

function findTraj(x,r_moon,v_moon,idx0,idxf, uu0 = [0.01*rand(6) for k = 1:idxf-idx0-1])
    N = idxf - idx0
    n = 12 # number of states
    m = 6 # number of controls

    #Model contains any paramters needed
    model  = R3BP{Float64}(earth.μ, moon.μ)

    ## EDIT THESE TO CHANGE SENSITIVITY TO INITIAL STATE OR CONTROL
    #Weights
    Q = 1.0*Diagonal(I,n) #Cost of intermediate states (make large if you want the trajectory to follow really closely)
    Qf = 1.0*Diagonal(I,n) #Cost of final state
    R = 1.0*Diagonal(I,m) #Cost of control (make large if you don't mind if the trajectory deviates as long as the control is small)

    #costfuns setting the cost of varying from the initial guess trajectory at each time step
    costfuns = map(1:N) do k #Setting
        xf = x[:,idx0+k]
        rfm = r_moon[:,idx0+k]
        vfm = v_moon[:,idx0+k]
        LQRCost(Q,R,[xf;rfm;vfm])
    end

    #Objective function made of the sum of all costfuns
    obj = Objective(costfuns)

    # #constraints = Constraints(N) # define constraints at each time step
    # for k = 1:N-1
    #     constraints[k] += bnd
    # end
    # constraints[N] += goal

    x0 = [x[:,idx0+1]; r_moon[:,idx0+1]; v_moon[:,idx0+1]] #initial state
    xf = [x[:,idxf]; r_moon[:,idxf]; v_moon[:,idxf]] #final state
    t0 = t[idx0+1] #initial time
    tf = t[idxf] #final time
    prob = Problem(model, obj, xf, tf, x0=x0, t0=t0) # construct problem
    initial_controls!(prob,uu0) # initialize problem with controls

    # solver = solve!(prob, ALTROSolverOptions{Float64}())
    solver = iLQRSolver(prob)

    TrajectoryOptimization.solve!(solver)

    xx = states(solver)
    uu = controls(solver)

    return xx, uu
end

function optim(xx0, et0, etf, uu0 = [0.01*rand(3) for k = 1:70740])
    M,N = size(xx0)
    n = 6 # number of states
    m = 3 # number of controls

    # Set up interpolations of lunar position and attitude
    dt = 600 #sec
    range = et0:dt:etf+dt
    r_moon_coarse = [spkpos("moon",et,"J2000","none","earth")[1] for et = range]
    r_moon = scale(interpolate(r_moon_coarse, BSpline(Linear())),range)
    R_coarse = [pxform("MOON_PA", "J2000", et) for et = range]
    R = scale(interpolate(R_coarse, BSpline(Linear())),range)

    #Model contains any paramters needed
    model  = EPHEM{Float64}(earth.μ, moon.μ, sun.μ, r_moon, R)

    #Weights
    Q = 1.0e-2*Diagonal(I,n) #Cost of intermediate states
    Qf = 1.0e-2*Diagonal(I,n) #Cost of final state
    R = 1.0e2*Diagonal(I,m) #Cost of control

    #costfuns setting the cost of varying from the initial guess trajectory at each time step
    costfuns = map(1:N) do k
        xf = xx0[:,k]
        LQRCost(Q,R,xf)
    end

    #Objective function made of the sum of all costfuns
    obj = Objective(costfuns)

    # #constraints = Constraints(N) # define constraints at each time step
    # for k = 1:N-1
        # constraints[k] += bnd
    # end
    # constraints[N] += goal

    x0 = xx0[:,1] #initial state
    xf = xx0[:,end] #final state
    prob = Problem(model, obj, xf, etf, x0=x0, t0=et0) # construct problem
    initial_controls!(prob,uu0) # initialize problem with controls

    solver = iLQRSolver(prob)

    TrajectoryOptimization.solve!(solver)
    xx = states(solver)
    uu = controls(solver)
    return xx,uu
end

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
