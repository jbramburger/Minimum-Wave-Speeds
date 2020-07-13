%--------------------------------------------------------------------------
% Application of the surface method for upper bounds on the minimum wave 
% speed to the modified Fisher equation:
%
%      u_t = u_xx + u^m(1-u),    m = 1,2,3,...
%
% This code is associated to the paper "Minimum wave speeds in monostable 
% reaction-diffusion equations: sharp bounds by polynomial optimization" by 
% Jason J. Bramburger and David Goluskin (2020). This script is used to 
% create the data in Tables 1 and 2. 
%--------------------------------------------------------------------------

% Access YALMIP and mosek directories
addpath(genpath('YALMIP-master')) 
addpath(genpath('mosek'))

% Clean workspace
clear all
close all
clc

format long

% Differential equation parameter
m = 2; %f(u) = (1-u)*u^m

% Bounding method parameter
d = 1; %Degree of N(u)

%Bisection Method
cleft = 0;
cright = 2;
while abs(cright - cleft) >= 1e-5

    cmid = 0.5*(cleft + cright);
    flag = surface(cmid,m,d);
    
    if flag == 0
       cright = cmid; 
    else
       cleft = cmid;
    end
end

%Print Results
fprintf('The upper bound on the minimum speed for m = %d is %f found using degree %d polynomials.\n',m,cmid,d)

%% 
function flag = surface(c,m,d)

% SDP variable declarations
sdpvar u
z = sdpvar(2,1);

% Auxiliary function
[N, cN] = polynomial(u,d);

% S Procedure
d2 = d+m;
[s11, s11c] = polynomial(u,d2,0);
[s12, s12c] = polynomial(u,d2,0);
[s22, s22c] = polynomial(u,d2,0);
S = [s11, s12; s12, s22]; %Symmetric 2x2 matrix of auxillary functions

% Function defintions and integral terms
r = 1/(m+2)*u^(m+2) - 1/(m+1)*u^(m+1);
intN = int(N,u,0,u);
Q = [-c*intN+r, N; N, 2];

% Constraints
cons = [sos(z'*(Q - u*(1-u)*S)*z), sos(z'*S*z)];

% SOS Solver
ops = sdpsettings('solver','mosek','verbose',0);
sol = solvesos(cons,[],ops,[cN; s11c; s12c; s22c]);

% Return whether solvesos failed or succeeded
flag = sol.problem;

end
























