function [rv_J2000,t] = CR3BP_2_J2000(rv_CR3BP,theta,pos_sec,vel_sec,mu,origin)
% rv_J2000 = CR3BP_2_J2000(rv_CR3BP,t,origin)
% converts state vector from rotating frame in CR3BP to inertial frame,
% cenetered at origin in J2000
% Inputs:
%     rv_rot = [6x1] state vector [rx;ry;rz;vx;vy;vz] in CR3BP rotating frame
%     theta = (scalar) radians rotated
%     origin = (string) 'BARY', 'PRIM', 'SEC' origin of inertial frame
% Outputs
%     rv_inert = [6x1] state vector [rx;ry;rz;vx;vy;vz] in inertial frame

[n,m] = size(rv_CR3BP);
if n ~= 6
    if m == 6
        rv_CR3BP = rv_CR3BP';
        m = n;
    else
        error("rv_rot should be 6xN")
    end
end

global RUNIT VUNIT TUNIT


if nargin < 6;      origin = "PRIM";    end
if nargin < 5;      global mu;          end
if isempty(mu);     mu = 3.0404233e-06; end %Value for Europa/Jupiter System


if length(theta) == 1
    theta = theta*ones(1,m);
elseif iscolumn(theta)
    theta = theta';
end

if upper(origin) == "BARY"
elseif upper(origin) == "PRIM"
    rv_CR3BP = rv_CR3BP + [mu;0;0;0;0;0];
elseif upper(origin) == "SEC"
    rv_CR3BP = rv_CR3BP - [1-mu;0;0;0;0;0];
else
    error('origin must be BARY, PRIM, or SEC')
end

rv_CR3BP(1:3,:) = rv_CR3BP(1:3,:)*RUNIT;
rv_CR3BP(4:6,:) = rv_CR3BP(4:6,:)*VUNIT;

rv_J2000 = zeros(size(rv_CR3BP));
for i = 1:length(theta)
    R = pos_sec(:,i);
    V = vel_sec(:,i);
    x_hat = R/norm(R);
    z_hat = cross(R,V)/norm(cross(R,V));
    y_hat = cross(z_hat,x_hat);
    theta_dot = norm(cross(R,V))/norm(R)^2;
    B = [theta_dot*y_hat, -theta_dot*x_hat, zeros(3,1)];
    C = [x_hat, y_hat, z_hat];
    
    A = [C zeros(3,3);B,C];
    rv_J2000(:,i) = A*rv_CR3BP(:,i);
end
   
t = theta*TUNIT;

end