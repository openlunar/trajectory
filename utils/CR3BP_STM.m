function [statedot] = CR3BP_STM(t, state, mu)
%Equations of motion for the Circular Restricted Three-Body Problem

if nargin < 3; global mu; end

if length(state) == 42
    phi = state(1:36);
    PHI = reshape(phi, 6, 6);
    rv = state(37:42);
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
    
    omgxx= 1+(mu2/r5)*(3*(x+mu)^2)+(mu/R5)*(3*(x-mu2)^2)-(mu2/r3+mu/R3);
    omgyy= 1+(mu2/r5)*(3* y^2    )+(mu/R5)*(3* y^2     )-(mu2/r3+mu/R3);
    omgzz=   (mu2/r5)*(3* z^2    )+(mu/R5)*(3* z^2     )-(mu2/r3+mu/R3);
    
    omgxy= 3*y*  (mu2*(x+mu)/r5+mu*(x-mu2)/R5);
    omgxz= 3*z*  (mu2*(x+mu)/r5+mu*(x-mu2)/R5);
    omgyz= 3*y*z*(mu2       /r5+mu        /R5);
    
    F     =[   0     0     0     1     0	 0 ;
        0     0     0     0 	   1 	 0 ;
        0	 0     0     0     0     1 ;
        omgxx omgxy omgxz     0     2 	 0 ;
        omgxy omgyy omgyz    -2     0 	 0 ;
        omgxz omgyz omgzz     0	   0	 0 ];
    
    
    mu1 = mu2;
    X = rv;
    F(4,1)      = 1 - mu1*(1/r13 - 3*xx1^2/r15) - mu*(1/r23 - 3*xx2^2/r25);
    F(4,2)      = 3*mu1*xx1*X(2)/r15 + 3*mu*xx2*X(2)/r25;
    F(4,3)      = 3*mu1*xx1*X(3)/r15 + 3*mu*xx2*X(3)/r25;
    F(5,1)      = F(4,2);
    F(5,2)      = 1 - mu1*(1/r13 - 3*X(2)^2/r15) - mu*(1/r23 - 3*X(2)^2/r25);
    F(5,3)      = 3*mu1*X(2)*X(3)/r15 + 3*mu*X(2)*X(3)/r25;
    F(6,1)      = F(4,3);
    F(6,2)      = F(5,3);
    F(6,3)      = - mu1*(1/r13 - 3*X(3)^2/r15) - mu*(1/r23 - 3*X(3)^2/r25);
    
    PHIdot = F * PHI;
    phidot = reshape(PHIdot, 36, 1);
    
    statedot        = zeros(42,1);
    statedot(1:36)  = phidot;
    statedot(37)    = vx;
    statedot(38)    = vy;
    statedot(39)    = vz;
    statedot(40)    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    statedot(41)    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
    statedot(42)    =  -(mu2* z    /r3) -(mu* z     /R3);

elseif length(state) == 20
    phi = state(1:16);
    PHI = reshape(phi, 4, 4);
    rv = state(17:20);
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
    
    omgxx= 1+(mu2/r5)*(3*(x+mu)^2)+(mu/R5)*(3*(x-mu2)^2)-(mu2/r3+mu/R3);
    omgyy= 1+(mu2/r5)*(3* y^2    )+(mu/R5)*(3* y^2     )-(mu2/r3+mu/R3);
    
    omgxy= 3*y*  (mu2*(x+mu)/r5+mu*(x-mu2)/R5);

    
    F =[0       0       1       0;
        0       0       0       1;
        omgxx   omgxy   0       2;
        omgxy   omgyy   -2      0];
    
    
    mu1 = mu2;
    X = rv;
    F(3,1)      = 1 - mu1*(1/r13 - 3*xx1^2/r15) - mu*(1/r23 - 3*xx2^2/r25);
    F(3,2)      = 3*mu1*xx1*X(2)/r15 + 3*mu*xx2*X(2)/r25;
    F(4,1)      = F(3,2);
    F(4,2)      = 1 - mu1*(1/r13 - 3*X(2)^2/r15) - mu*(1/r23 - 3*X(2)^2/r25);
    
    PHIdot = F * PHI;
    phidot = reshape(PHIdot, 16, 1);
    
    statedot        = zeros(20,1);
    statedot(1:16)  = phidot;
    statedot(17)    = vx;
    statedot(18)    = vy;
    statedot(19)    = x-(mu2*(x+mu)/r3) -(mu*(x-mu2)/R3) + 2*vy;
    statedot(20)    = y-(mu2* y    /r3) -(mu* y     /R3) - 2*vx;
end

% % rvdot = [rv(4);
% %     rv(5);
% %     rv(6);
% %     (-mu/(r^3))*x;
% %     (-mu/(r^3))*y;
% %     (-mu/(r^3))*z];
% % 
% % G = [(-mu/(r^3)),0,0;
% %     0,(-mu/(r^3)),0;
% %     0,0,(-mu/(r^3))];
% % 
% % %     G = (mu/(r^5))*[3*x^2-r^2, 3*x*y,3*x*z;...
% % %         3*y*x, 3*y^2-r^2, 3*y*z;...
% % %         3*z*x, 3*z*y, 3*z^2-r^2];
% % 
% % H = zeros(3);
% % 
% % F = [zeros(3), eye(3);
% %     G, H];
% % 
% % PHIdot_mat = F*PHI_mat;
% % PHIdot = [reshape(PHIdot_mat,36,1);
% %     rvdot];

end

