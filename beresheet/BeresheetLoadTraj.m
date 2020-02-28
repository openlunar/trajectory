clc
close all
clear all
%%
addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
addpath('~/trajectory/beresheet/')
addpath('~/StanfordMATLAB/')

%%
setEarthMoonGlobal

%% Read in data from .txt file
xx = readtable("BeresheetTrajectory.txt");

%% Extract dates and trajectory data
dates = table2array(xx(:,1));
infmt = 'dd MMM yyyy HH:mm:ss.SSS';
datetimes = datetime(dates,'InputFormat',infmt);
t0 = datenum(datetimes(1));
t = datenum(datetimes) - t0;
t = t*24*60*60; %[sec]

x = table2array(xx(:,2:end));

%% rotating back into the xy plane
Rx = rotx(-23.5);
Ry = roty(-5);
plot_rv([Ry*Rx*x(:,1:3)';x(:,4:6)'])

%%
x(:,1:3) = x(:,1:3)./RUNIT;
x(:,4:6) = x(:,4:6)./VUNIT;
t = t/TUNIT;
xrot = inert2rot(x',t)';
xt = [xrot,t];

plot_CR3BP
plot_rv(xrot)

%% finding the maneuvers and discontinuities

r   = zeros(1,length(x));
v   = zeros(1,length(x));
dv  = zeros(1,length(x));
for i = 1:(length(x))
    r(i) = norm(x(i,1:3));
    v(i) = norm(x(i,4:6));
    if i > 1
        dv(i) = norm(x(i,4:6) - x(i-1,4:6));
    end
end

figure()
plot(r)
xlabel('#')
ylabel('r [km]')
grid on
title('Distance from Earth')

figure()
plot(v)
xlabel('#')
ylabel('v [km/s]')
grid on
title('Velocity Relative to Earth')

figure()
hold on
plot(r/max(r))
plot(v/max(v))
plot(dv)
grid on
legend('r','v','dv')
title('Normalized r and v with dv showing discontinuities')

%% Finding number of maneuvers
EarthPeri = islocalmin(r(1:58000)); %before we get to the Moon
sum(EarthPeri)

[dv_sort, dv_ind] = sort(dv);
discontinuities = dv_ind(end-1:end)
%% Finding periapses and apoapses
% It looks like there are three maneuvers to get to the Moon
% So there are four sections of the orbit

RE = 6.378136e3; %[km]
mu = 3.986004328969392e5; %[km^3/s^2] for Earth

rA = zeros(1,4); hA = rA; %[km]
rP = zeros(1,4); hP = rP; %[km]
vA = zeros(1,4); %[km/s]
vP = zeros(1,4); %[km/s]
e = zeros(1,4);
T = zeros(1,4);

rA(1) = max(r(1:4000)); hA(1) = rA(1) - RE; %[km]
rP(1) = min(r(1:4000)); hP(1) = rP(1) - RE; %[km]
vA(1) = min(v(1:4000)); %[km/s]
vP(1) = max(v(1:4000)); %[km/s]
e(1) = (rA(1) - rP(1))/(rA(1) + rP(1));
a(1) = rA(1)/(1+e(1));
T(1) = 2*pi*sqrt(a(1)^3/mu)

rA(2) = max(r(12000:18000)); hA(2) = rA(2) - RE; %[km]
rP(2) = min(r(12000:18000)); hP(2) = rP(2) - RE; %[km]
vA(2) = min(v(12000:18000)); %[km/s]
vP(2) = max(v(12000:18000)); %[km/s]
e(2) = (rA(2) - rP(2))/(rA(2) + rP(2));
a(1) = rA(1)/(1+e(1));

rA(3) = max(r(20000:30000)); hA(3) = rA(3) - RE; %[km]
rP(3) = min(r(20000:30000)); hP(3) = rP(3) - RE; %[km]
vA(3) = min(v(20000:30000)); %[km/s]
vP(3) = max(v(20000:30000)); %[km/s]
e(3) = (rA(3) - rP(3))/(rA(3) + rP(3));

rA(4) = max(r(40000:54000)); hA(4) = rA(4) - RE; %[km]
rP(4) = min(r(40000:54000)); hP(4) = rP(4) - RE; %[km]
vA(4) = min(v(40000:54000)); %[km/s]
vP(4) = max(v(40000:54000)); %[km/s]
e(4) = (rA(4) - rP(4))/(rA(4) + rP(4));




%%
% setearth
% plot_prim('b')
% hold on
%%
figure()
comet3(x(:,1)/RUNIT, x(:,2)/RUNIT, x(:,3)/RUNIT)