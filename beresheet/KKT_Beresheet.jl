using ForwardDiff, Plots, SparseArrays
using LinearAlgebra, MAT, TrajectoryOptimization

function CR3BPdynamics!(X) #Three body dynamics in Earth/Moon System
    m = 7.34767309e22 # [kg] Mass of the moon
    M = 5.9722e24 # [kg] Mass of the earth
    mu = m/(M+m)

    x = 1*X[1];
    y = 1*X[2];
    z = 1*X[3];
    vx = 1*X[4];
    vy = 1*X[5];
    vz = 1*X[6];

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

    Ẋ = zeros(eltype(X),6)
    # @show typeof(Ẋ)
    Ẋ[1]    = 1*vx;
    Ẋ[2]    = 1*vy;
    Ẋ[3]    = 1*vz;
    Ẋ[4]    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    Ẋ[5]    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
    Ẋ[6]    =  -(mu2* z    /r3) -(mu* z     /R3);
    return Ẋ
end

function f(x̄) # fourth order Runge-Kutta Expansion #NOTE: I had to put everything into the x̄ vector so I can easily take derivatives
    if length(x̄) < 6
        error("f(x̄): x̄ should be at least 6x1")
    elseif length(x̄) == 6
        x, u, h = x̄, zeros(6), 0.01 #default time step
        println("default time step")
    elseif length(x̄) == 7
        x, u, h = x̄[1:6], zeros(6), x̄[7]
    elseif length(x̄) == 8
        error("f(x̄): x̄ should not be 8x1")
    elseif length(x̄) == 9
        x, u, h = x̄[1:6], [zeros(3);x̄[7:9]], 0.01 #default time step
        println("default time step")
    elseif length(x̄) == 10
        x, u, h = x̄[1:6], [zeros(3);x̄[7:9]], x̄[10]
    else
        error("f(x̄): x̄ should be smaller than 10x1")
    end
    k1 = CR3BPdynamics!(x) + u
    k2 = CR3BPdynamics!(x + h/2*k1) + u
    k3 = CR3BPdynamics!(x + h/2*k2) + u
    k4 = CR3BPdynamics!(x + h*k3) + u
    return x + h*(1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4)
end

Af = x̄ -> ForwardDiff.jacobian(f, x̄); #Jacobian of Linearlize Dynamics
J(x,u) = (x-x₀)'*(x-x₀) + ones(6)'*[I(3); -I(3)]*u #Objective function

function c(X̄)
    N = Int64((length(X̄)+4)/10)
    # N = 2

    x̄ = X̄[1:6N]
    ū = X̄[6N+1:9N-3]
    dt̄ = X̄[9N-2:end]
    # c = Vector{Float64}(undef, 6*(N-1))
    c = Array{Any,1}(undef,6*(N-1))
    # c = zeros(6*(N-1))
    for n = 1:N-1
        c[n*6-5:n*6] = x̄[(n+1)*6-5:(n+1)*6] - f([x̄[n*6-5:n*6]; ū[n*3-2:n*3]; dt̄[n]])
    end
    return c
end

# function c(xₙ₊₁, xₙ, uₙ, dtₙ)
#     x = x̄[1:6]
#     u = x̄[7:9]
#     dt = x̄[10]
#     # c = Vector{Float64}(undef, 6*(N-1))
#
#     c = x[(n+1)*6-5:(n+1)*6] - f([x[n*6-5:n*6]; u[n*3-2:n*3]; dt[n]])
#
# end

Ac = x̄ -> ForwardDiff.jacobian(c, x̄); #Jacobian of Linearlize Dynamics

# function BeresheetKKT()
    # Load in the starting trajectory
cd("/home/rexlab/trajectory/beresheet/")
file = matopen("BeresheetTrajectory.mat")
xt = read(file,"BeresheetTraj")
close(file)
N,M = size(xt)
x̄₀ = reshape(xt[:,1:6]',N*6,1)
t̄₀ = xt[:,7];
dt̄₀ = t̄₀[2:end] - t̄₀[1:end-1]
ū₀ = zeros(3*(N-1))

#Regularization into CR3BP coords
RUNIT = 384400 #[km] Distance from Earth to Moon in CR3BP
TUNIT = 27.32*24*60*60/2π # [sec] Sidereal Period of the moon is about 27.32 days
VUNIT = RUNIT/TUNIT
for i = 1:N
    x̄₀[6i-5:6i-3] /= RUNIT
    x̄₀[6i-2:6i] /= VUNIT
end
t̄₀ /= TUNIT
dt̄₀ /= TUNIT


# model = Model(CR3BPdynamics!,6,1)
# model_d = rk4(model)
# end

function nonLinearKKT(x̄₀; b=1, iter=20)
    N = length(x̄₀)
    ū = zeros(3*(N-1))


    H1(x̄) =[I Ac(x̄)]
    H2(x̄) =[Ac(x̄)' 0]
    H̄(x̄) = [H1(x̄);H2(x̄)]

    x̄ = x̄₀
    # println("x = ", x)
    # println("|x| = ", norm(x))

    # # plot(x)
    # x₁ =  Array{Float64,1}(undef, iter)
    # x₂ =  Array{Float64,1}(undef, iter)
    # x₁[1] = x[1]
    # x₂[1] = x[2]

    # z₁ = zeros(iter)
    # z₂ = zeros(iter)

    δx = [1,1]

    # while norm(δx) > 0.001
    for i = 1:iter-1
        # global x
        # global x₁
        # global x₂

        Ā = H̄(x)

        b̄ = [-x + x₀;-c(x)]

        δxλ = Ā\b̄
        δx = δxλ[1:2]
        λ = δxλ[3]

        α = 1
        newpoint = false
        while newpoint == false
            if (Φ(x+α*δx,x₀) <= Φ(x,x₀))
                x = x + α*δx
                newpoint = true
            elseif α == 1
                δx̂ = -J(x)*inv(J(x)'*J(x))*c(x + δx)
                if (Φ(x+δx+δx̂,x₀) <= Φ(x,x₀))
                    x = x + δx + δx̂
                    newpoint = true

                    z = x + δx
                    z₁[i] = z[1]
                    z₂[i] = z[2]
                    println("I used the quadratic step")
                else
                    α = α/2
                end
            else
                α = α/2
            end
        end

        # x = x + α*δx
        x₁[i+1] = x[1]
        x₂[i+1] = x[2]
        println("x = ", x)
        println("|x| = ", norm(x))
        # plot!(x)

    end
end

#Using TrajectoryOptimization Package
