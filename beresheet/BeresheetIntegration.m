clc
clear all
close all

addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')

%%
setEarthMoonGlobal
load BeresheetTrajectory.mat

%%
tf = (BeresheetTraj(end,end) - BeresheetTraj(36630,end))*24*60*60/TUNIT;
x0 = BeresheetTraj(38000,1:6);
x0(1:3) = x0(1:3)/RUNIT;
x0(4:6) = x0(4:6)/(VUNIT);
x0 = inert2rot(x0);
[tt, xx] = ode78e(@(t,y) CR3BP(t,y),0,100*tf,x0,1e-12);

%%
xx_inert = rot2inert(xx,tt);
plot_rv(xx_inert)

%%
plot_CR3BP
plot_rv(xx)