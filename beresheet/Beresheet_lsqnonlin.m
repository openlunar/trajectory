clc
clear all
close all

addpath('~/trajectory/beresheet/')
addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
% load BeresheetCR3BP1000.mat
setEarthMoonGlobal

load BeresheetCR3BP.mat
%%
% x0 = reshape(xCR3BP_1000,6000,1);
N = 50;
% for i = 1
i=7;
x0 = reshape(rv_CR3BP(:,N*(i-1)+1:N*i),6*N,1);

% t0 = 

u0 = zeros(3*(N-1),1);
J = @(x) (x(1:6*N) - x0)'*(x(1:6*N) - x0) + x(6*N+1:9*N-3)'*x(6*N+1:9*N-3); %l2 norm square just for simplicity right now

%%



% PROBLEM.objective = J;
% PROBLEM.x0 = [x0;u0];
% PROBLEM.nonlcon = C;
% PROBLEM.solver = 'fmincon';
% X = fmincon(PROBLEM);
OPTIONS = optimoptions('fmincon', 'MaxFunctionEvaluations',1e4);
x = [x0;u0];
t = theta(1:N)';

X = lsqnonlin(@(x) NONLCON2(x),x);

%%
xx = reshape(X(1:6*N),6,N);
uu = reshape(X(6*N+1:9*N-3),3,N-1);

figure(1)
plot_rv(xx)
plot_rv(rv_CR3BP(:,N*(i-1)+1:N*i),'r')
figure(2)

subplot(4,1,1)
plot(uu(1,:))
title('u_x')
subplot(4,1,2)
plot(uu(2,:))
title('u_y')
subplot(4,1,3)
plot(uu(3,:))
title('u_z')
subplot(4,1,4)
for j = 1:N-1
    unorm(j) = norm(uu(:,j));
end
plot(unorm)
title('||u||')
% end 
%%
% [C,Ceq] = NONLCON(x)

%%
% pwd
% cd StanfordMATLAB
open AA279_HW2_2.m

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

function dX = trapezoid(X)
N = (length(X)+3)/9;
xk = reshape(X(1:(N-1)*6),6,N-1);
xk1 = reshape(X(7:N*6),6,N-1);
uk = reshape(X(N*6+1:N*9-6),3,N-2);
uk1 = reshape(X(N*6+4:N*9-3),3,N-2);
for i = 1:N-1
    h = 60/3.7570e+05;
    if i == N-1
        fk = CR3BP(0,xk(:,i)) + [0;0;0;uk1(:,i-1)];
        fk1 = CR3BP(0,xk1(:,i));
    else
        fk = CR3BP(0,xk(:,i)) + [0;0;0;uk(:,i)];
        fk1 = CR3BP(0,xk1(:,i)) + [0;0;0;uk1(:,i)];
    end
    dX(:,i) = h/2*(fk+fk1);
end
dX = reshape(dX,(N-1)*6,1);
end

function xk1 = hermite(xk, hk)
f = @(tk,xk,uk) CR3BP(tk,xk) + [0;0;0;uk];

% xc = 1/2*(xk+xk1) + hk/8*(fk-fk1)
% xk1 = xk + hk/6*( + 4*(xc) + )

N = (length(X)+3)/9;
xk = reshape(X(1:(N-1)*6),6,N-1);
xk1 = reshape(X(7:N*6),6,N-1);
uk = reshape(X(N*6+1:N*9-6),3,N-2);
uk1 = reshape(X(N*6+4:N*9-3),3,N-2);
for i = 1:N-1
    h = 60/3.7570e+05;
    if i == N-1
         fk = CR3BP(0,xk(:,i)) + [0;0;0;uk1(:,i-1)];
        fk1 = CR3BP(0,xk1(:,i));
    else
        fk = CR3BP(0,xk(:,i)) + [0;0;0;uk(:,i)];
        fk1 = CR3BP(0,xk1(:,i)) + [0;0;0;uk1(:,i)];
    end
    dX(:,i) = h/2*(fk+fk1);
end
dX = reshape(dX,(N-1)*6,1);
end


function [C, Ceq] = NONLCON(x)
% Ceq = x(7:6000) - RK4(x);
N = (length(x)+3)/9;
Ceq = x(7:6*N) - x(1:6*(N-1)) - trapezoid(x);
C = [];
end

function Ceq = NONLCON2(x)
% Ceq = x(7:6000) - RK4(x);
N = (length(x)+3)/9;
Ceq = x(7:6*N) - x(1:6*(N-1)) - trapezoid(x);
end