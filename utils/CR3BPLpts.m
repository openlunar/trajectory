function Req = CR3BPLpts(tol)
% Req = CR3BPLpts(tol)
% 
% computes the 5 Lagrange points in the CR3BP.
% 
% inputs:
% tol       :tolerance to use for iterative solver (difference between two
%            successive approximations)
% 
% outputs:
% Req       :array of 5 position vectors for L1-L5
% 
% Record of Revision
% Date          Programmer          Description of Changes
% 11/06/2012    Brian D. Anderson   Original Code
% 
% Brian D. Anderson
% bdanders@usc.edu
% University of Southern California
% Los Angeles, CA


% access global variables
global mu

% set default tolerance
if nargin<1
    tol         = 1e-15;
end

% solve for L1
alpha           = (mu/3*(1 - mu))^(1/3); %initial guess
dalpha          = 1;

% iteratively solve for alpha
while abs(dalpha)>tol
    alphaprev   = alpha;
    alpha       = (mu*(1 - alpha)^2/...
        (3 - 2*mu - alpha*(3 - mu - alpha)))^(1/3);
    dalpha      = alpha - alphaprev;
end

% compute x coordinate for L1
L1x             = 1 - mu - alpha;

% solve for L2
beta            = (mu/3*(1 - mu))^(1/3); %initial guess
dbeta           = 1;

% iteratively solve for beta
while abs(dbeta)>tol
    betaprev    = beta;
    beta        = (mu*(1 + beta)^2/...
        (3 - 2*mu + beta*(3 - mu + beta)))^(1/3);
    dbeta       = beta - betaprev;
end

% compute x coordinate for L2
L2x             = 1 - mu + beta;

% solve for L3
gamma           = -7/12*mu + 1; %initial guess
dgamma          = 1;

% iteratively solve for gamma
while abs(dgamma)>tol
    gammaprev   = gamma;
    gamma       = ((1- mu)*(1 + gamma)^2/...
        (1 + 2*mu + gamma*(2 + mu + gamma)))^(1/3);
    dgamma      = gamma - gammaprev;
end

% compute x coordinate for L3
L3x             = - mu - gamma;

% determine triangular points:
L4x             = 0.5 - mu;
L4y             = sqrt(3)/2;
L5x             = L4x;
L5y             = - L4y;

Req             = [L1x 0 0;L2x 0 0 ;L3x 0 0 ;L4x L4y 0;L5x L5y 0]';