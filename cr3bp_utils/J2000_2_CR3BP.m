function [rv_CR3BP,theta] = J2000_2_CR3BP(rv_J2000,t,pos_sec,vel_sec,mu,origin)
% rv_J2000 = CR3BP_2_J2000(rv_CR3BP,t,origin)
% converts state vector from rotating frame in CR3BP to inertial frame,
% cenetered at origin in J2000
% Inputs:
%     rv_J2000 = [6x1] state vector [rx;ry;rz;vx;vy;vz] in CR3BP rotating frame
%     theta = (scalar) radians rotated
%     origin = (string) 'BARY', 'PRIM', 'SEC' origin of inertial frame
% Outputs
%     rv_inert = [6x1] state vector [rx;ry;rz;vx;vy;vz] in inertial frame

[n,m] = size(rv_J2000);
if n ~= 6
    if m == 6
        rv_J2000 = rv_J2000';
        m = n;
    else
        error("rv_rot should be 6xN")
    end
end

if nargin < 6;      origin = "PRIM";    end
if nargin < 5;      global mu;          end
if isempty(mu);     mu = 3.0404233e-06; end %Value for Europa/Jupiter System

global RUNIT VUNIT TUNIT

if length(t) == 1
    t = t*ones(1,m);
elseif iscolumn(t)
    t = t';
end


rv_CR3BP = zeros(size(rv_J2000));
for i = 1:length(t)
    R = pos_sec(:,i);
    V = vel_sec(:,i);
    x_hat = R/norm(R);
    z_hat = cross(R,V)/norm(cross(R,V));
    y_hat = cross(z_hat,x_hat);
    theta_dot = norm(cross(R,V))/norm(R)^2;
    B = [theta_dot*y_hat, -theta_dot*x_hat, zeros(3,1)];
    C = [x_hat, y_hat, z_hat];
    
    J2000_RUNIT = norm(R);
    J2000_TUNIT = 1/sqrt(J2000_RUNIT^3/mu);
    J2000_VUNIT = J2000_RUNIT/J2000_TUNIT;
    
    A = [C zeros(3,3);B,C];
    rv_CR3BP(:,i) = inv(A)*rv_J2000(:,i);
    rv_CR3BP(1:3,i) = rv_CR3BP(1:3,i)/J2000_RUNIT;
    rv_CR3BP(4:6,i) = rv_CR3BP(4:6,i)/J2000_VUNIT;
end


if upper(origin) == "BARY"
elseif upper(origin) == "PRIM"
    rv_CR3BP = rv_CR3BP - [mu;0;0;0;0;0];
elseif upper(origin) == "SEC"
    rv_CR3BP = rv_CR3BP + [1-mu;0;0;0;0;0];
else
    error('origin must be BARY, PRIM, or SEC')
end

theta = t/TUNIT;
    
   
end