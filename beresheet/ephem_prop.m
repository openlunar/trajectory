function rvdot = ephem_prop(t, rv, mu)
%Ephemeris Propagator
load('BeresheetMoonSunEphem.mat')

if nargin < 3; global mu; end

if length(rv) == 6
    x = rv(1);
    y = rv(2);
    z = rv(3);
    vx = rv(4);
    vy = rv(5);
    vz = rv(6);
    
    mu2 = 1-mu;
    
    % r1 = (( x+mu )^2 + y^2 + z^2)^(1/2); %distance to m1, LARGER MASS
    % r2 = ((x-1+mu)^2 + y^2 + z^2)^(1/2); %distance to m2, smaller mass
    
    r2= (x+mu )^2 + y^2 + z^2;	% r: distance to m1, LARGER MASS
    R2= (x-mu2)^2 + y^2 + z^2;	% R: distance to m2, smaller mass
    r3= r2^1.5; r5= r2^2.5;
    R3= R2^1.5; R5= R2^2.5;
    
    xx1 = x + mu;
    xx2 = x-mu2;
    r12 = r2;
    r22 = R2;
    r13 = r3;
    r23 = R3;
    r15 = r5;
    r25 = R5;
    
    mu1 = mu2;
    
    rvdot        = zeros(6,1);
    rvdot(1)    = vx;
    rvdot(2)    = vy;
    rvdot(3)    = vz;
    rvdot(4)    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    rvdot(5)    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
    rvdot(6)    =  -(mu2* z    /r3) -(mu* z     /R3);
    
elseif length(rv) == 4
    x = rv(1);
    y = rv(2);
    vx = rv(3);
    vy = rv(4);
    
    mu2 = 1-mu;
    
    % r1 = (( x+mu )^2 + y^2 + z^2)^(1/2); %distance to m1, LARGER MASS
    % r2 = ((x-1+mu)^2 + y^2 + z^2)^(1/2); %distance to m2, smaller mass
    
    r2= (x+mu )^2 + y^2;	% r: distance to m1, LARGER MASS
    R2= (x-mu2)^2 + y^2;	% R: distance to m2, smaller mass
    r3= r2^1.5; r5= r2^2.5;
    R3= R2^1.5; R5= R2^2.5;
    
    xx1 = x + mu;
    xx2 = x-mu2;
    r12 = r2;
    r22 = R2;
    r13 = r3;
    r23 = R3;
    r15 = r5;
    r25 = R5;
    
    mu1 = mu2;
    
    rvdot       = zeros(4,1);
    rvdot(1)    = vx;
    rvdot(2)    = vy;
    rvdot(3)    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    rvdot(4)    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
else
    error('Wrong size state, should be 6x1 or 4x1')
end

%Changelog
%Date           Programmer              Action
%02/20/2020     Jared T. Blanchard      File created
