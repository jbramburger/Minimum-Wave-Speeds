%--------------------------------------------------------------------------
% Application of the volume method for lower bounds on the minimum wave 
% speed to the modified Fisher equation:
%
%      u_t = u_xx + u^m(1-u),    m = 1,2,3,...
%
% This code is associated to the paper "Minimum wave speeds in monostable 
% reaction-diffusion equations: sharp bounds by polynomial optimization" by 
% Jason J. Bramburger and David Goluskin (2020). This script is used to 
% create the data in Figure 2 and Table 2. 
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

% Bounding method parameters
d = 3; %Degree of V(u,v)
lambda = 10;

%Bisection Method
cleft = 0;
cright = 2;
while abs(cright - cleft) >= 1e-5

    cmid = 0.5*(cleft + cright);
    flag = volume(cmid,m,d,lambda);
    
    if flag == 0
       cleft = cmid; 
    else
       cright = cmid;
    end
end

%Print Results
fprintf('The lower bound on the minimum speed for m = %d is %f found using degree %d polynomials.\n',m,cmid,d)

%%
function flag = volume(c,m,d,lambda)

% SDP variables
sdpvar u v

% Epsilon value
eps = 1e-4;

% Auxiliary function
[V, cV] = polynomial([u v], d);

% S Procedure
d2 = d+m;
[s1, c1] = polynomial([u v], d2);
[s2, c2] = polynomial([u v], d2);
[s3, c3] = polynomial(u, d2);

% Derivatives
dVdu = jacobian(V,u);
dVdv = jacobian(V,v);

% Replacements
Vv0 = replace(V,v,0);

%Constraints
cons = [];
cons = [cons, replace(V,[u v], [0 -eps]) == 0, replace(V, [u v], [1 0]) == 0];
cons = [cons, sos(lambda*(dVdu*v - c*dVdv*v - dVdv*(1-u)*u^m) - V - u*(1-u)*s1 + v*s2)]; 
cons = [cons, sos(-Vv0 - eps*(1-u) - u*(1-u)*s3)];
cons = [cons, sos(s1), sos(s2), sos(s3)];

%SOS Solver
ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
sol = solvesos(cons,[],ops,[cV;c1; c2; c3]);

%Return whether solvesos failed or succeeded
flag = sol.problem;

end














