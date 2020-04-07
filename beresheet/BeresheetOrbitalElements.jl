using SatelliteDynamics
include("/home/rexlab/trajectory/utils/OrbitDynamics.jl")

load_gravity_model("/home/rexlab/trajectory/kernels/STU_MoonTopo720.gfc")

rv2oe()
