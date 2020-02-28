clear all
close all
clc

%%
prev_dir = pwd;
cd('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/jpl/spice')
miceAddPath

addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
addpath('~/trajectory/beresheet/')
addpath('~/StanfordMATLAB/')
setEarthMoonGlobal

%% Read in data from .txt file
xx = readtable("BeresheetTrajectory.txt");

%% Extract dates and trajectory data
dates = table2array(xx(:,1));
infmt = 'dd MMM yyyy HH:mm:ss.SSS';
datetimes = datetime(dates,'InputFormat',infmt);
t0 = datenum(datetimes(1));
t = datenum(datetimes) - t0; %starting at zero
t = t*24*60*60; %[sec]

x = table2array(xx(:,2:end))';


%%

% Construct a meta kernel, "standard.tm”, which will be used to load the needed
% generic kernels: "naif0011.tls," "de421.bsp,” and "pck00010.tpc.”
% Load the generic kernels using the meta kernel, and a Cassini spk.
cspice_furnsh( { 'kernels/standard.tm', 'kernels/030201AP_SK_SM546_T45.bsp'} )
% Define the number of divisions of the time interval.
STEP = length(x);
et = cspice_str2et( {char(datetimes(1)), char(datetimes(end))} );
%      date = 'Thu Mar 20 12:53:29 PST 1997';
%
%      %
%      % Convert a string to ephemeris time (ET).
%      %
%      et = cspice_str2et( date );

times = (0:STEP-1) * ( et(2) - et(1) )/STEP + et(1); %seconds from Jan 1, 2000, 12 hr

[pos_moon,~]= cspice_spkpos( 'MOON', times, 'J2000', 'NONE', 'EARTH' ); %position in km
vel_moon = zeros(size(pos_moon));
vel_moon(:,1) = (pos_moon(:,2) - pos_moon(:,1))/(times(2) - times(1));
for i = 2:length(pos_moon)-1
    vel_moon(:,i) = (pos_moon(:,i+1) - pos_moon(:,i-1))/(times(i+1) - times(i-1));
end
vel_moon(:,end) = (pos_moon(:,end) - pos_moon(:,end-1))/(times(end) - times(end-1));


[pos_sun,~]= cspice_spkpos( 'SUN', times, 'J2000', 'NONE', 'EARTH' );
vel_sun = zeros(size(pos_sun));
vel_sun(:,1) = (pos_sun(:,2) - pos_sun(:,1))/(times(2) - times(1));
for i = 2:length(pos_sun)-1
    vel_sun(:,i) = (pos_sun(:,i+1) - pos_sun(:,i-1))/(times(i+1) - times(i-1));
end
vel_sun(:,end) = (pos_sun(:,end) - pos_sun(:,end-1))/(times(end) - times(end-1));


cspice_kclear

cd(prev_dir)

%%
save('~/trajectory/beresheet/BeresheetMoonSunEphem.mat','x','pos_moon','pos_sun','vel_moon','vel_sun','times')

%% Plot the resulting trajectory with ephemeris motion
figure()
hold on
plot_rv([pos_moon;pos_moon],'b')
% plot_rv([pos_sun;pos_sun],'y')
plot_rv(x,'r')
for i = 1:100:length(x)
    h1 = plot3(pos_moon(1,i),pos_moon(2,i),pos_moon(3,i),'bo');
%     h2 = plot3(pos_sun(1,i), pos_sun(2,i), pos_sun(3,i),'yo');
    h3 = plot3(x(1,i), x(2,i), x(3,i),'ro');
    drawnow
    delete(h1)
%     delete(h2)
    delete(h3)
end

%% Transform to CR3BP
[rv_CR3BP,theta] = J2000_2_CR3BP(x,t,pos_moon,vel_moon);
figure()
plot_CR3BP
plot_rv(rv_CR3BP)

%% Compare the TUNIT used in J2000_2_CR3BP with CR3BP TUNIT
figure()
plot(t'./theta)
hold on
plot(TUNIT*ones(size(t)))
title('TUNIT')
legend('Ephem','CR3BP')

%%
save('~/trajectory/beresheet/BeresheetCR3BP.mat','rv_CR3BP','theta','pos_moon','vel_moon')
%% Try integrating a portion to see where it goes
% x0 = rv_CR3BP(:,38660);
x0 = rv_CR3BP(:,1);
[tt_integrated, rv_integrated] = ode78e(@(t,y) CR3BP(t,y),0,1,x0,1e-12);
plot_CR3BP
plot_rv(rv_CR3BP,'k')
plot_rv(rv_integrated,'r')
plot3(x0(1),x0(2),x0(3),'rx')

%% convert back to J2000
[rv_J2000, t_J2000] = CR3BP_2_J2000(rv_CR3BP,theta,pos_moon,vel_moon);
figure()
plot_rv(rv_J2000,'r')
plot_rv(x,'k')
legend('reconverted','original')

%% compare velocities
figure()
plot_rv(rv_J2000(4:6,:),'r')
plot_rv(x(),'k')
legend('reconverted','original')



%% look closer at the moon
pos_rel_moon = pos_moon - x(:,1:3)';

figure(2)
hold on
plot_rv([pos_rel_moon;pos_rel_moon])
plot_sphere(MOON.radius,[0,0,0],2,[.6,.6,.6])
axis([-1e4,1e4,-1e4,1e4,-1e4,1e4])

% plot_rv(x,'r')
% for i = 1:100:length(x)
%     h1 = plot3(pos_moon(1,i),pos_moon(2,i),pos_moon(3,i),'bo');
%     h2 = plot3(pos_sun(1,i), pos_sun(2,i), pos_sun(3,i),'yo');
%     h3 = plot3(x(i,1), x(i,2), x(i,3),'ro');
%     drawnow
%     delete(h1)
%     delete(h2)
%     delete(h3)
% end



for i = 1:length(pos_rel_moon)
a(i) = norm(pos_rel_moon(:,i)) - MOON.radius;
end
figure()
plot(a)

%%

r = x(38996,1:3) %38996 is the one when the sun is closest to the vernal equinox
v = x(38996,1:3)

%58877 is the one when the moon is closest to the vernal equinox
r = x(38668,1:3)
v = x(38668,1:3)

%cd back into the previous directory

