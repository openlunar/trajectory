clc
clear all
close all

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
X = fmincon(@(x) J(x),x,[],[],[],[],[],[],@(x) NONLCON(x))

%%
xx = reshape(X(1:6000),6,1000);
uu = reshape(X(6001:9000),3,1000);

plot_rv(xx)

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

function [C, Ceq] = NONLCON(x)
Ceq = x(7:6000) - RK4(x);
C = [];
end
