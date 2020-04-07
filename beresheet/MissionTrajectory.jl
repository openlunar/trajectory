dir = "/home/rexlab/trajectory/"
include(string(dir,"utils/OrbitDynamics.jl"))
include(string(dir,"utils/OptimizationSetup.jl"))

# Setup Spice Kernel
furnsh(string(dir,"kernels/naif0012.tls"))
furnsh(string(dir,"kernels/de432s.bsp"))
furnsh(string(dir,"kernels/moon_pa_de421_1900-2050.bpc"))
furnsh(string(dir,"kernels/moon_080317.tf"))

et0 = utc2et("2019-02-22T02:18:40")
etf = utc2et("2019-04-12T05:18:00")

# R = pxform("MOON_PA", "J2000", et0)
# S = sxform("MOON_PA", "J2000", et0)
rm = spkpos("moon",et0,"J2000","none","earth")[1]
# spkezr("moon",et0,"J2000","none","earth")

file = matopen(string(dir,"beresheet/fullBeresheetTraj.mat"))
xx0 = read(file,"x")
uu0 = read(file,"u")
et = read(file,"et")
t = read(file,"t")
close(file)

plot(xx0[1,:],xx0[2,:],title="Trajectory - Actual Launch Date")
uu0[:,5600:6200] = zeros(3,601)
uu0[:,36400:36860] = zeros(3,461)
uu0_norm = [norm(uu0[:,i]) for i = 1:70740]
plot(uu0_norm,title="Maneuvers - Actual Launch Date")

idx = 70000
r = xx0[1:3,idx]
rm = spkpos("moon",et0+(idx-1)*60,"J2000","none","earth")[1]
R = pxform("MOON_PA", "J2000", et0+(idx-1)*60)
a1 = -moon.μ*(r-rm)/norm(r-rm)^3
a2 = -accel_gravity(1000*(r-rm),R,10,10)

#Load in Beresheet trajectory
month = 27.321661*86400
furnsh(string(dir,"beresheet/fullBeresheetTraj.bsp"))
et0b = utc2et("2019-02-22T02:20:48.183")
etfb = utc2et("2019-04-12T05:18:05.187")
boddef("beresheet",-5441)
xx_beresheet = [spkezr("beresheet",et,"J2000","none","earth")[1] for et = et0b:60:etfb]
xx_beresheet = convert2array(xx_beresheet)



# Change the date here
xx,uu = optim(xx_beresheet,et0+month,etf+month)



xx_array = convert2array(xx)
plot(xx_array[1,:],xx_array[2,:],title="Trajectory - 1 Month Later",label="")
uu_norm = [norm(uu[i]) for i = 1:length(uu)]
uu_norm[5600:6200] = zeros(601)
uu_norm[36400:36860] = zeros(461)
plot(uu_norm,title="Maneuvers - 1 Month Later",label="u")

furnsh(string(dir,"beresheet/mission.bsp"))
et1 = utc2et("2022-06-16T21:45:52.556")
et2 = utc2et("2022-06-21T15:55:52.556")
r1 = spkezr("-5440",et1,"J2000","none","earth")[1]
r2 = spkezr("-5440",et2,"J2000","none","earth")[1]

xx_mission = [spkezr("-5440",et,"J2000","none","earth")[1] for et = et1:60:et2]

xx_mission_array = convert2array(xx_mission)
plot(xx_mission_array[1,:],xx_mission_array[2,:])
