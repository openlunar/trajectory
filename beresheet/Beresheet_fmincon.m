clc
clear all
close all

addpath('~/trajectory/beresheet/')
addpath('~/jared711@stanford.edu/STANFORD/Research/OrbitalResearch/misc/util/')
% load BeresheetCR3BP1000.mat
setEarthMoonGlobal

load BeresheetCR3BP.mat
load BeresheetTrajectory.mat
rv_J2000 = BeresheetTraj(:,1:6)';
%%
% x0 = reshape(xCR3BP_1000,6000,1);
N = 200;
% for i = 1
i = 7;
% x0 = reshape(rv_CR3BP(:,N*(i-1)+1:N*i),6*N,1);
x0 = reshape(rv_J2000(:,N*(i-1)+1:N*i),6*N,1);

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
% X = fmincon(@(x) J(x),x,[],[],[],[],[],[],@NONLCON,OPTIONS);
X = lsqnonlin(@(x) NONLCON3(x),x);

%%
xx = reshape(X(1:6*N),6,N);
uu = reshape(X(6*N+1:9*N-3),3,N-1);

figure(1)
plot_rv(xx)
plot_rv(rv_J2000(:,N*(i-1)+1:N*i),'r')
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

function dX = hermite(X)
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
        
        xc = 1/2*(xk(:,i)+xk1(:,i)) + h/8*(fk-fk1);
        fc = CR3BP(0,xc) + [0;0;0;uk1(:,i-1)/2];
    else
        fk = CR3BP(0,xk(:,i)) + [0;0;0;uk(:,i)];
        fk1 = CR3BP(0,xk1(:,i)) + [0;0;0;uk1(:,i)];
        
        xc = 1/2*(xk(:,i)+xk1(:,i)) + h/8*(fk-fk1);
        fc = CR3BP(0,xc) + [0;0;0;uk(:,i)+uk1(:,i)/2];
    end
    dX(:,i) = h/6*(fk+4*fc+fk1);
end
dX = reshape(dX,(N-1)*6,1);
end

function dX = hermite2(X)
N = (length(X)+3)/9;
xk = reshape(X(1:(N-1)*6),6,N-1);
xk1 = reshape(X(7:N*6),6,N-1);
uk = reshape(X(N*6+1:N*9-6),3,N-2);
uk1 = reshape(X(N*6+4:N*9-3),3,N-2);
for i = 1:N-1
    h = 60/3.7570e+05;
    if i == N-1
        fk = EarthGrav(0,xk(:,i)) + [0;0;0;uk1(:,i-1)];
        fk1 = EarthGrav(0,xk1(:,i));
        
        xc = 1/2*(xk(:,i)+xk1(:,i)) + h/8*(fk-fk1);
        fc = EarthGrav(0,xc) + [0;0;0;uk1(:,i-1)/2];
    else
        fk = EarthGrav(0,xk(:,i)) + [0;0;0;uk(:,i)];
        fk1 = EarthGrav(0,xk1(:,i)) + [0;0;0;uk1(:,i)];
        
        xc = 1/2*(xk(:,i)+xk1(:,i)) + h/8*(fk-fk1);
        fc = EarthGrav(0,xc) + [0;0;0;uk(:,i)+uk1(:,i)/2];
    end
    dX(:,i) = h/6*(fk+4*fc+fk1);
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
% Ceq = x(7:6*N) - x(1:6*(N-1)) - trapezoid(x);
Ceq = x(7:6*N) - x(1:6*(N-1)) - hermite(x);
end

function Ceq = NONLCON3(x)
% Ceq = x(7:6000) - RK4(x);
N = (length(x)+3)/9;
% Ceq = x(7:6*N) - x(1:6*(N-1)) - trapezoid(x);
Ceq = x(7:6*N) - x(1:6*(N-1)) - hermite2(x);
end

function rvdot = EarthGrav(t,rv)
%Gravity from Earth only
r = rv(1:3);
mu = 3.98600751696e+05; %km^3/s^2
rvdot = zeros(6,1);
rvdot(1:3) = rv(4:6);
rvdot(4:6) = -mu*r/norm(r)^3;
end
