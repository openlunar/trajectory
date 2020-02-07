clc
clear all
close all

addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
load BeresheetCR3BP1000.mat
setEarthMoonGlobal

%%
x0 = reshape(xCR3BP_1000,6000,1);
u0 = zeros(3000,1);
J = @(x) (x(1:6000) - x0)'*(x(1:6000) - x0) + x(6001:9000)'*x(6001:9000); %l2 norm square just for simplicity right now

%%



% PROBLEM.objective = J;
% PROBLEM.x0 = [x0;u0];
% PROBLEM.nonlcon = C;
% PROBLEM.solver = 'fmincon';
% X = fmincon(PROBLEM)
x = [x0;u0];
X = fmincon(@(x) J(x),x,[],[],[],[],[],[],@NONLCON);

%%
xx = reshape(X(1:6000),6,1000);
uu = reshape(X(6001:9000),3,1000);

plot_rv(xx)
%%

diff = max(max(xx - reshape(x0,6,1000)))

%%
trapezoid(x)

%%
% [C,Ceq] = NONLCON(x)

%%
function x1 = RK4(xu)
global mu
x0 = reshape(xu(1:5994),6,999);
u0 = reshape(xu(6001:8997),3,999);
h = (4*pi + pi)/1000;
for i = 1:999
    k1 = CR3BP(0,x0(:,i),mu) + [0;0;0;u0(:,i)];
    k2 = CR3BP(0,x0(:,i)+h*k1/2,mu) + [0;0;0;u0(:,i)];
    k3 = CR3BP(0,x0(:,i)+h*k2/2,mu) + [0;0;0;u0(:,i)];
    k4 = CR3BP(0,x0(:,i)+h*k3,mu) + [0;0;0;u0(:,i)];
    x1(:,i) = x0(:,i) + h*(k1/6 + k2/3 + k3/3 + k4/6);
end
x1 = reshape(x1,5994,1);
end

function xk1 = trapezoid(X)
global mu
xk = reshape(X(1:5994),6,999);
xk1 = reshape(X(7:6000),6,999);
uk = reshape(X(6001:8997),3,999);
uk1 = reshape(X(6004:9000),3,999);
h = (4*pi + pi)/1000;
for i = 1:999
    fk = CR3BP(0,xk(:,i),mu) + [0;0;0;uk(:,i)];
    fk1 = CR3BP(0,xk1(:,i),mu) + [0;0;0;uk1(:,i)];
    
    xk1(:,i) = xk(:,i) + h/2*(fk+fk1);
end
xk1 = reshape(xk1,5994,1);
end

function [C, Ceq] = NONLCON(x)
% Ceq = x(7:6000) - RK4(x);
Ceq = x(7:6000) - trapezoid(x);
C = [];
end
