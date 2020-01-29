using ForwardDiff, Plots, SparseArrays
using LinearAlgebra, MAT, TrajectoryOptimization

function cubicSpline(x₀,t₀,k)
#x(t) = at^3 + bt^2 + ct + d where a,b,c,d ∈ ℜ⁶
    M,N = size(x₀) #let's assume that x ∈ ℜᴹˣᴺ
                    # M is the size of the state
                    # N is the number of data points
    # T₁ = ceil(N/k)
    # T₂ = floor(N/k)
    #
    # if N != N%k * T₁ + (k-N%k)*T₂
    #     error("won't have enough knot points")
    # end
    SSᵣ = zeros(M,1);

    q = 0
    t = t₀
    abcd = zeros(4M,k)
    x = zeros(M,k)
    for i = 1:k
        if i <= N%k
            T = Int64(ceil(N/k))
        else
            T = Int64(floor(N/k))
        end

        A = zeros(M*T,M*4)
        b = zeros(M*T,1)
        j₀ = 0
        for j = 1:T
            q += 1
            A[j₀+1:j₀+M,:] = [I(6)*t[q]^3 I(6)*t[q]^2 I(6)*t[q] I(6)*1]
            b[j₀+1:j₀+M] = x₀[:,q]
            j₀ += M
        end
        abcd[:,i] = A\b
        a = abcd[1:M,i]
        b = abcd[M+1:2M,i]
        c = abcd[2M+1:3M,i]
        d = abcd[3M+1:4M,i]
        x[:,i] = a*t[q]^3 + b*t[q]^2 + c*t[q] + d

        SSᵣ += sqrt.((x₀[:,q] - x[:,i]).^2)
    end

    return x, SSᵣ
end

cd("/home/rexlab/trajectory/beresheet/")
file = matopen("BeresheetTrajectory.mat")
xt = read(file,"BeresheetTraj")
close(file)
N,M = size(xt)
x̄₀ = reshape(xt[:,1:6]',6,70742)
t̄₀ = xt[:,7];
dt̄₀ = t̄₀[2:end] - t̄₀[1:end-1]
ū₀ = zeros(3*(N-1))

x, SSᵣ = cubicSpline(x̄₀,t̄₀,1000)
plot(x̄₀[1,:],x̄₀[2,:])
plot(x[1,:],x[2,:])

Juno.@enter cubicSpline(x̄₀,t̄₀,1000)
