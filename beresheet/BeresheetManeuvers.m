clc
clear all
close all

% I need both the util and the StanfordMATLAB directories on my path

setearth
PRIM = EARTH;
RUNIT = EARTH.sm;

r = EARTH.radius

%% Setting up the planned manuevers as shown in the Beresheet Paper 

% SEP.StartTime = 
SEP.Apogee = 70000 + r;
SEP.Perigee = 215 + r;

% AM1.StartTime = %gotta remember the conversion to UTCG
AM1.Apogee = 68715 + r; %[km]
AM1.Perigee = 215 + r; %[km]

% AM2.StartTime = %gotta remember the conversion to UTCG
AM2.Apogee = 68710 + r; %[km]
AM2.Perigee = 600 + r; %[km]

% PM1.StartTime = %gotta remember the conversion to UTCG
PM1.Apogee = 117300 + r; %[km]
PM1.Perigee = 530 + r; %[km]

% PM2.StartTime = %gotta remember the conversion to UTCG
PM2.Apogee = 275000 + r; %[km]
PM2.Perigee = 1600 + r; %[km]

% PM3.StartTime = %gotta remember the conversion to UTCG
PM3.Apogee = 390000 + r; %[km]
PM3.Perigee = 1500 + r; %[km]

% OPM.StartTime = %gotta remember the conversion to UTCG
OPM.Apogee = 392000 + r; %[km]
OPM.Perigee = 1800 + r; %[km]

%% First let's try to plot the Separation orbit

SEP.a = (SEP.Apogee + SEP.Perigee)/2
SEP.e = ap2e(SEP.Apogee, SEP.Perigee)
SEP.nu = 0;
SEP.omega = 0; %[deg]
SEP.i = 0; %[deg]
SEP.OMEGA = 0; %[deg]

[r_ECI, v_ECI] = oe2rv(SEP.a*1000,SEP.e,SEP.nu,...
    SEP.omega,SEP.i,SEP.OMEGA);

y0 = [r_ECI/1000*RUNIT;v_ECI/1000*VUNIT];

[tt, xx] = ode78e(@(t,y) CR3BPbrian(t,y),0,10,y0,1e-12,0);

plot_prim('b')
plot_rv(xx(1:5,:))

