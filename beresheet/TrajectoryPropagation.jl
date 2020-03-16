using DifferentialEquations
include("/home/rexlab/trajectory/utils/OrbitDynamics.jl")
include("/home/rexlab/trajectory/utils/OrbitDynamics.jl")

file = matopen("/home/rexlab/trajectory/beresheet/BeresheetMoonSunEphem.mat")
r_moon = read(file,"pos_moon")
v_moon = read(file,"vel_moon")
r_sun = read(file,"pos_sun")
v_sun = read(file,"vel_sun")
x = read(file,"x")
t = read(file,"times")
close(file)

rv, theta = J2000_2_CR3BP(x,t,[r_moon;v_moon])

rv0 = rv[:,100]
tspan = (theta[1],10*theta[1000])
μ = 0.012150577032698
prob = ODEProblem(CR3BPdynamics!,rv0,tspan,μ)
alg = TsitPap8()
sol = solve(prob,alg)
sol2 = solve(prob)

plot_earth(CR3BP=true)
p=plot(sol,
    vars=(1,2,3),
    xlabel = "X [NON]",
    ylabel = "Y [NON]",
    zlabel = "Z [NON]",
    xlims = (-.5,.5),
    ylims = (-.5,.5),
    zlims = (-.5,.5),
    aspect_ratio = :equal
)
plot!(sol2,
    vars=(1,2,3),
    xlabel = "X [NON]",
    ylabel = "Y [NON]",
    zlabel = "Z [NON]",
    xlims = (-.5,.5),
    ylims = (-.5,.5),
    zlims = (-.5,.5),
    aspect_ratio = :equal
)

plot(sol,vars=(1,2))


anim = animate(sol,"/home/rexlab/trajectory/beresheet/trajprop.gif")


## Try using ode8()
