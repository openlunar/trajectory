include(string(pwd(),"/utils/OrbitDynamics.jl"))
include(string(pwd(),"/utils/OptimizationSetup.jl"))


# et0 = utc2et("2019-02-22T02:18:40") #Beresheet initial launch date
# etf = utc2et("2019-04-12T05:18:00") #Beresheet landing date

# Example code showing how to call Spice functions
# R = pxform("MOON_PA", "J2000", et0)
# S = sxform("MOON_PA", "J2000", et0)
# rm = spkpos("moon",et0,"J2000","none","earth")[1]
# rv = spkezr("moon",et0,"J2000","none","earth")

#Load in initial guess from a mat file
# file = matopen(string(pwd(),"/beresheet/fullBeresheetTraj.mat"))
# xx0 = read(file,"x") #state
# uu0 = read(file,"u") #control policy
# et = read(file,"et") #ephemeris times
# t = read(file,"t") #time in seconds starting from 0
# close(file)

# #Plotting the loaded Trajectory
# plot(xx0[1,:],xx0[2,:],title="Trajectory - Actual Launch Date")
# uu0[:,5600:6200] = zeros(3,601)
# uu0[:,36400:36860] = zeros(3,461)
# uu0_norm = [norm(uu0[:,i]) for i = 1:70740]
# plot(uu0_norm,title="Maneuvers - Actual Launch Date")

# #Sanity checking the lunar gravity model acceleration
# idx = 70000
# r = xx0[1:3,idx]
# rm = spkpos("moon",et0+(idx-1)*60,"J2000","none","earth")[1]
# R = pxform("MOON_PA", "J2000", et0+(idx-1)*60)
# a1 = -moon.μ*(r-rm)/norm(r-rm)^3
# a2 = -accel_gravity(1000*(r-rm),R,10,10)

# Try using the Beresheet initial conditions, but a month later
month = 27.321661*86400 # one month in [sec]
furnsh(string(pwd(),"/beresheet/fullBeresheetTraj.bsp")) #Load in Beresheet trajectory
et0b = utc2et("2019-02-22T02:20:48.183") #Actual Beresheet launch date
etfb = utc2et("2019-04-12T05:18:05.187") #Actual Beresheet landing date
boddef("beresheet",-5441) #define ephemeris object -5441 as Beresheet
xx_beresheet = [spkezr("beresheet",et,"J2000","none","earth")[1] for et = et0b:60:etfb]
xx_beresheet = convert2array(xx_beresheet) #array needs to be fed into optim()



# Change the date here
xx,uu = optim(xx_beresheet,et0+month,etf+month) #just one month in the future
# xx,uu = optim(xx_beresheet,et0+12*month,etf+12*month) #one year in the future

# Plotting the results
xx_array = convert2array(xx) #Easier to plot if in array form
plot(xx_array[1,:],xx_array[2,:],title="Trajectory - 1 Month Later",label="")
uu_norm = [norm(uu[i]) for i = 1:length(uu)]
uu_norm[5600:6200] = zeros(601)
uu_norm[36400:36860] = zeros(461)
plot(uu_norm,title="Maneuvers - 1 Month Later",label="u")

# #Setup the result of patched_conic.py in the optimizer
# furnsh(string(pwd(),"/beresheet/mission.bsp"))
# et1 = utc2et("2022-06-16T21:45:52.556")
# et2 = utc2et("2022-06-21T15:55:52.556")
# r1 = spkezr("-5440",et1,"J2000","none","earth")[1]
# r2 = spkezr("-5440",et2,"J2000","none","earth")[1]
#
# xx_mission = [spkezr("-5440",et,"J2000","none","earth")[1] for et = et1:60:et2]
#
# xx_mission_array = convert2array(xx_mission)
# plot(xx_mission_array[1,:],xx_mission_array[2,:])
