clc
clear all
close all

addpath('jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
setEarthMoonGlobal

t0 = 0;
tf = 10;

x0 = [0.5;0;0;0;0.5;0];
PHI0 = eye(6);
X0 = [reshape(PHI0,36,1);x0];

[tt,xx] = ode78e(@(t,y) CR3BP(t,y),t0,tf,x0);
figure()
plot_CR3BP
plot_rv(xx)
x = xx(end,:)

[ttSTM,xxSTM] = ode78e(@(t,y) CR3BP_STM(t,y),t0,tf,X0);
figure()
plot_CR3BP
plot_rv(xxSTM)
xSTM = xxSTM(end,37:42)
%%
PHI = reshape(xxSTM(end,1:36),6,6)

dx0 = [1e-3;0;0;0;0;0];
[tt1,xx1] = ode78e(@(t,y) CR3BP(t,y),t0,tf,x0 + dx0);
x = xx(end,:)
x1 = xx1(end,:)
dx = xx1(end,:) - xx(end,:)
error = (PHI*dx0)' - dx
