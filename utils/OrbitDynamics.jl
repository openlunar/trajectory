using Plots, Colors, LinearAlgebra
plotlyjs()

function xyz_extract(rv)
    n,m = size(rv)
    if n == 6 || n == 7
        x,y,z,vx,vy,vz = rv[1:6,:]
        return x,y,z,vx,vy,vz
    elseif n == 4
        x,y,vx,vy = rv[1:4,:]
        return x,y,vx,vy
    elseif m == 6 || m == 7
        rv = rv'
        x,y,z,vx,vy,vz = rv[1:6,:]
        return x,y,z,vx,vy,vz
    elseif m == 4
        rv = rv'
        x,y,vx,vy = rv[1:4,:]
        return x,y,vx,vy
    end
    return r,v
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

function rv_extract(rv)
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
    end
    return r,v
end

function R2BPdynamics(t,rv,mu=398600)
    r,v = rv_extract(rv)
    u = zeros(size(v))
    rdot = v
    vdot = mu/norm(r)^3 * r + u
    return [rdot;vdot]
end


function CR3BPdynamics(t,rv,mu = 0.012150577032698) #Three body dynamics in Earth/Moon System
    x,y,z,vx,vy,vz = xyz_extract(rv)

    mu2 = 1-mu;
    r1 = (( x+mu )^2 + y^2 + z^2)^(1/2); #distance to m1, LARGER MASS
    r2 = ((x-1+mu)^2 + y^2 + z^2)^(1/2); #distance to m2, smaller mass

    r2= (x+mu )^2 + y^2 + z^2; # r: distance to m1, LARGER MASS
    R2= (x-mu2)^2 + y^2 + z^2; # R: distance to m2, smaller mass
    r3= r2^1.5; r5= r2^2.5;
    R3= R2^1.5; R5= R2^2.5;

    xx1 = x + mu;
    xx2 = x-mu2;
    r12 = r2;
    r22 = R2;
    r13 = r3;
    r23 = R3;
    r15 = r5;
    r25 = R5;

    mu1 = mu2;

    rvdot = zeros(eltype(rv),6)
    # @show typeof(Ẋ)
    rvdot[1:3] = [vx,vy,vz]
    rvdot[4]    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    rvdot[5]    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
    rvdot[6]    =  -(mu2* z    /r3) -(mu* z     /R3);
    return rvdot
end

function plot_rv(rv)
    r,v = rv_extract(rv)
    p = plot!(r[1,:],r[3,:],r[3,:])
    return p
end

function plot_sphere(r=1,c=[0,0,0],col='r',n=100)
    u = range(0,stop=2π,length=n)
    v = range(0,stop=π,length=n)
    x = cos.(u) * sin.(v)';
    y = sin.(u) * sin.(v)';
    z = ones(n) * cos.(v)';
    const_color = cgrad( [ RGB{Float64}(0.,0.,1.) for _ in 1:2 ] )
    s = surface!(r*x,r*y,r*z,colorbar=false,color=const_color,width=1)
    return s
end
