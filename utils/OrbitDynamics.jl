using Plots
plotlyjs()
using Colors
using LinearAlgebra
using MAT
using Interpolations
using StaticArrays
using SPICE
using SatelliteDynamics
load_gravity_model(string(dir,"kernels/STU_MoonTopo720.gfc")) #load in lunar gravity field rather than Earth's
using TrajectoryOptimization
import TrajectoryOptimization.discrete_dynamics


##  Contains various useful functions for orbital dynamics, including
#   dynamics for R2BP and CR3BP as well as plotting and frame conversion

mutable struct body
    name
    color
    μ # gravitational parameter  [km³/s²]
    m # mass                     [kg]
    a # semimajor axis           [km]
    ω # rotation rate            [rad/s]
    T # period                   [sec]
    r # radius                   [km]
    e # eccentricity             [NON]
end

earth = body("Earth",[0.,0.,1.],398600,5.97237e24,1.495978740473789e8 ,1.160576151685878e-5,3.155814910224000e7, 6371.0084,0)
moon = body("Moon",[.1,.1,.1],4902.801,7.34767309e22,384400,2.661666501955582e-6,2.3606208e6,1737.5,0.0554)
sun = body("Sun",[0.6,0.6,0.],1.32712440042e11, 1.988500e30, 0, 0, 0, 6.957e5, 0)

mutable struct system
    μ
end

earth_moon = system(0.012150577032698)

function xyz_extract(rv)
    n,m = try n,m = size(rv)
    catch
        size(rv)[1],1
    end
    if n == 7
        x,y,z,vx,vy,vz,t = rv[1:7,:]
        return x,y,z,vx,vy,vz,t
    elseif n == 6
        x,y,z,vx,vy,vz = rv[1:6,:]
        return x,y,z,vx,vy,vz
    elseif n == 4
        x,y,vx,vy = rv[1:4,:]
        return x,y,vx,vy
    elseif n == 3
        x,y,z = rv[1:3,:]
        return x,y,z
    elseif m == 7
        rv = rv'
        x,y,z,vx,vy,vz,t = rv[1:7,:]
        return x,y,z,vx,vy,vz,t
    elseif m == 6
        rv = rv'
        x,y,z,vx,vy,vz = rv[1:6,:]
        return x,y,z,vx,vy,vz
    elseif m == 4
        rv = rv'
        x,y,vx,vy = rv[1:4,:]
        return x,y,vx,vy
    elseif m == 3
        rv = rv'
        x,y,z = rv[1:3,:]
        return x,y,z
    else
        error("xyz_extract(rv): rv has incorrect dimensions")
    end
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

function semimajor(rA,rP)
    a = (rA + rP)/2
    return a
end

function eccentricity(rA,rP)
    e = (rA-rP)/(rA+rP)
end

function vel(a,r,μ=398600)
    v = sqrt(2*μ/r - μ/a)
    return v
end

function apoPeriVel(rA,rP)
    a = semimajor(rA,rP)
    vA = vel(a,rA)
    vP = vel(a,rP)
    return vA,vP
end

function rv2oe(rv,μ = 3.986e5)
    # Outputs orbital elements for a given state [r;v]
    r = rv[1:3]
    v = rv[4:6]
    # specific angular momentum
    h = cross(r,v)
    # directions of angular momentum
    W = h/norm(h)

    i = atan(sqrt(W[1]^2+W[2]^2),W[3])

    if i == 0 || i == pi/2
        Ω = NaN
    else
        Ω = atan(W[1],-W[2])
    end

    p = norm(h)^2/μ
    a = (2/norm(r) - norm(v)^2/μ)^(-1)
    n = sqrt(μ/a^3)
    e = sqrt(1 - p/a)
    E = atan((r'*v)/(n*a^2),(1-norm(r)/a))
    ν = E2nu(E,e)
    u = atan(r[3]/sin(i),r[1]*cos(Ω) + r[2]*sin(Ω))
    ω = u-ν
    return a,e,i,Ω,ω,ν
end

function E2nu(E,e)
    #outputs true anomaly ν
    ν = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    while ν >= 2π
        ν = ν-2π;
    end
    while ν < 0
        ν = ν + 2π;
    end
    return ν
end

# Integrates the control policy into
function integ_control(uu,idx0=1,idxf=length(uu);dt=60)
    Δv = zeros(length(uu[1]))
    for n = idx0:idxf-1
        Δv .+= (uu[n+1] .+ uu[n])*dt/2
    end
    return Δv
end

function rv_extract!(rv)
    n,m = try n,m = size(rv)
    catch
        size(rv)[1],1
    end

    if n == 6 || n == 7
        r = rv[1:3,:]
        v = rv[4:6,:]
    elseif n == 4
        r = rv[1:2,:]
        v = rv[3:4,:]
    elseif m == 6 || m == 7
        rv = rv'
        r = rv[1:3,:]
        v = rv[4:6,:]
    elseif m == 4
        rv = rv'
        r = rv[1:2,:]
        v = rv[3:4,:]
    else
        error("rv_extract(rv): rv has incorrect dimensions")
    end
    return r,v
end

function R2BPdynamics!(rvdot,rv,μ,t)
    r,v = rv_extract!(rv)
    rvdot[1:3] = v
    rvdot[4:6] = μ/norm(r)^3 * r
end

function R3BPdynamics!(rvdot,rv,μ,t) # inplace dynamics
    r,v = rv_extract!(rv[1:6,:])
    rs,vs = rv_extract!(rv[7:12,:])
    rdot = v
    vdot = -μ[1]*r/norm(r)^3 - μ[2]*(r-rs)/norm(r-rs)^3
    rmdot = vm
    vmdot = -μ[1]*rm/norm(rm)^3
    rvdot = [rdot; vdot; rmdot; vmdot]
end

function R3BPdynamics(t,rv,μ=[earth.μ;moon.μ]) # inplace dynamics
    r,v = rv_extract!(rv[1:6,:])
    rm,vm = rv_extract!(rv[7:12,:])
    rdot = v
    vdot = -μ[1]*r/norm(r)^3 - μ[2]*(r-rm)/norm(r-rm)^3
    rmdot = vm
    vmdot = -μ[1]*rm/norm(rm)^3
    rvdot = [rdot; vdot; rmdot; vmdot]
end



function CR3BPdynamics!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = xyz_extract(rv)
    μ2 = 1-μ;
    r₁³= ( (x+μ )^2 + y^2 + z^2 )^1.5; # distance to m1, LARGER MASS
    r₂³= ( (x-μ2)^2 + y^2 + z^2 )^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4]   = x-(μ2*(x+μ)/r₁³) -(μ*(x-μ2)/r₂³) + 2*vy;
    rvdot[5]   = y-(μ2* y   /r₁³) -(μ* y    /r₂³) - 2*vx;
    rvdot[6]   =  -(μ2* z   /r₁³) -(μ* z    /r₂³);
end

function CR3BPdynamics(t,rv,μ=0.012150577032698) #Three body dynamics in Earth/Moon System
    rvdot = zeros(6,1)
    x,y,z,vx,vy,vz = xyz_extract(rv)
    μ2 = 1-μ;
    r₁³= ( (x+μ )^2 + y^2 + z^2 )^1.5; # distance to m1, LARGER MASS
    r₂³= ( (x-μ2)^2 + y^2 + z^2 )^1.5; # distance to m2, smaller mass
    rvdot[1:3] = [vx;vy;vz]
    rvdot[4]   = x-(μ2*(x+μ)/r₁³) -(μ*(x-μ2)/r₂³) + 2*vy;
    rvdot[5]   = y-(μ2* y   /r₁³) -(μ* y    /r₂³) - 2*vx;
    rvdot[6]   =  -(μ2* z   /r₁³) -(μ* z    /r₂³);
    return rvdot
end

function J2000_2_CR3BP(rv_J2000,t,rv_sec,μ=0.012150577032698)
    r,v = rv_extract!(rv_J2000)
    rs,vs = rv_extract!(rv_sec)
    rv_CR3BP = zeros(size(rv_J2000))
    theta = zeros(size(t))
    for i = 1:length(t)
        R = rs[:,i];
        V = vs[:,i];
        x_hat = R/norm(R);
        z_hat = cross(R,V)/norm(cross(R,V));
        y_hat = cross(z_hat,x_hat);
        theta_dot = norm(cross(R,V))/norm(R)^2;
        B = [theta_dot*y_hat -theta_dot*x_hat zeros(3)];
        C = [x_hat y_hat z_hat];

        J2000_RUNIT = norm(R);
        J2000_TUNIT = sqrt(J2000_RUNIT^3/398600);
        J2000_VUNIT = J2000_RUNIT/J2000_TUNIT;

        theta[i] = (t[i]-t[1])/J2000_TUNIT

        A = [C zeros(3,3);B C];
        rv_CR3BP[:,i] = inv(A)*rv_J2000[:,i];
        rv_CR3BP[1:3,i] = rv_CR3BP[1:3,i]/J2000_RUNIT;
        rv_CR3BP[4:6,i] = rv_CR3BP[4:6,i]/J2000_VUNIT;
    end
    rv_CR3BP = rv_CR3BP - [μ;0;0;0;0;0]*ones(1,length(t));
    return rv_CR3BP, theta
end

function plot_rv(rv;dim=3)
    r,v = rv_extract!(rv)
    if dim == 2
        p = plot!(r[1,:],r[2,:])
    elseif dim == 3
        p = plot!(r[1,:],r[2,:],r[3,:])
    else
        error("dim should be 2 or 3")
    end
    return p
end

function plot_sphere(r=1,c=[0,0,0],col='b',n=100)
    u = range(0,stop=2π,length=n)
    v = range(0,stop=π,length=n)
    x = cos.(u) * sin.(v)';
    y = sin.(u) * sin.(v)';
    z = ones(n) * cos.(v)';
    const_color = cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
    s = surface!(r*x,r*y,r*z,colorbar=false,color=const_color,width=1)
    return s
end

function plot_earth(;CR3BP=true,col='b',n=100)
    CR3BP ? r = earth.r/moon.a      : r = earth.r
    CR3BP ? c = [-earth_moon.μ,0,0] : c = [0,0,0]
    plot_sphere(r,c,col,n)
end

function plot_circle(r=1,c=[0,0,0],col='b',n=100)
    u = range(0,stop=2π,length=n)
    x = r*cos.(u)
    y = r*sin.(u)
    c = plot!(x,y)
    return c
end

function ode8(F,t0,tf,y0,tol=1e-6,etol=1e-13;h=(tf-t0)/1000)
    y = y0
    t = t0

    tout = t
    yout = y
    # τ = tol*max(norm(y,Inf),1)
    if tf > t0
        while (t < tf)
            if t+h > tf; h = tf-t; end
            #step forward in time
            t = t + h
            y = rk8(F,y,t,h)
            tout = [tout; t]
            yout = [yout y]
        end
    else
        h = -h
        while (t > tf)
            if t+h < tf; h = tf-t; end
            #step forward in time
            t = t + h
            y = rk8(F,y,t,h)
            tout = [tout; t]
            yout = [yout y]
        end
    end
    yout = [yout[:,i] for i = 1:length(tout)]
    return tout, yout

end



function rk8(F,y,t,h)
    α    = [ 2/27; 1/9; 1/6; 5/12; .5; 5/6; 1/6; 2/3; 1/3; 1; 0; 1 ];
    β    = [ 2/27       0       0      0        0         0       0         0     0      0     0 0 0;
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
                -1777/4100 0       0      -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1 0]'
    χ     = [0; 0; 0; 0; 0; 34/105; 9/35; 9/35; 9/280; 9/280; 0; 41/840; 41/840 ]; # Difference between two bottom layers of butcher tableau
    # ψ     = [1 0 0 0 0 0      0    0    0     0     1 -1     -1     ]';
    pow   = 1 / 8;

    f = y*zeros(1,13)

    f[:,1] = F(t,y)
    for j = 1:12
        f[:,j+1] = F(t + α[j]*h,y + h*f*β[:,j])
    end

    return y + h*f*χ
end

#HARGRAVES IMPLEMENTATION
function cubicSpline(x₀,t₀,k)
#x(t) = at^3 + bt^2 + ct + d where a,b,c,d ∈ ℜ⁶
    M,N = size(x₀) #let's assume that x ∈ ℜᴹˣᴺ
                    # M is the size of the state
                    # N is the number of data points
    # T₁ = ceil(N/k)
    # T₂ = floor(N/k)
    #``
    # if N != N%k * T₁ + (k-N%k)*T₂
    #     error("won't have enough knot points")
    # end
    SSᵣ = zeros(M,1);

    q = 1


    # k = knot*3 + 1
    numStages = Int64(floor((k-1)/3))

    abcd = zeros(4M,numStages)
    x = zeros(M,k)
    t = zeros(1,k)
    Nₛ = Array{Int64,1}(undef,numStages)
    Tₛ = zeros(1,numStages)
    for i = 1:numStages
        if i <= N%numStages
            Nₛ[i] = Int64(ceil(N/numStages))
        else
            Nₛ[i] = Int64(floor(N/numStages))
        end

        Tₛ[i] = t₀[Int64(q+Nₛ[i])] - t₀[q]

        τ = zeros(Int64(Nₛ[i]))
        for j = 1:Int64(Nₛ[i])
            τ[j] = (t₀[Int64(q+i)]-t₀[Int64(q+i-1)])/Tₛ[i]
        end

        A = zeros(M*Nₛ[i],M*4)
        b = zeros(M*Nₛ[i],1)
        j₀ = 0
        tₛ = 0
        for j = 1:Nₛ[i]
            q += 1
            A[j₀+1:j₀+M,:] = [I(6)*tₛ^3 I(6)*tₛ^2 I(6)*tₛ I(6)*1]
            b[j₀+1:j₀+M] = x₀[:,q]
            j₀ += M
            tₛ += τ[j]
        end
        abcd[:,i] = A\b
    end

    for i = 1:numStages
        a = abcd[1:M,i]
        b = abcd[M+1:2M,i]
        c = abcd[2M+1:3M,i]
        d = abcd[3M+1:4M,i]

        x[:,4*(i-1)+1] = a*(0*Tₛ[i]/3)^3 + b*(0*Tₛ[i]/3)^2 + c*(0*Tₛ[i]/3) + d
        x[:,4*(i-1)+2] = a*(1*Tₛ[i]/3)^3 + b*(1*Tₛ[i]/3)^2 + c*(1*Tₛ[i]/3) + d
        x[:,4*(i-1)+3] = a*(2*Tₛ[i]/3)^3 + b*(2*Tₛ[i]/3)^2 + c*(2*Tₛ[i]/3) + d

        # t[4*(i-1)+1] = T
        SSᵣ += sqrt.((x₀[:,q] - x[:,i]).^2)
    end
    return x, t, SSᵣ
end
