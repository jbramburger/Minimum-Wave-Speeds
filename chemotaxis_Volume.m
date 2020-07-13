%--------------------------------------------------------------------------
% Application of the volume method for lower bounds on the minimum wave 
% speed to the single-component reaction-diffusion equation:
%
%    u_t = (u^k u_x)_x - buu_x + u(1-u^q),   b > 0, k,q - positive integers    
%
% studied in the paper "Traveling wave solutions of a nonlinear reaction-
% diffusion-chemotaxis model for bacterial pattern formation" by M.B.A. 
% Mansour (2008)
%
% This code is associated to the paper "Minimum wave speeds in monostable 
% reaction-diffusion equations: sharp bounds by polynomial optimization" by 
% Jason J. Bramburger and David Goluskin (2020). This script is used to 
% create the data in Table 3. 
%--------------------------------------------------------------------------

% Access YALMIP and mosek directories
addpath(genpath('YALMIP-master')) 
addpath(genpath('mosek'))

% Clean workspace
clear all
close all
clc

format long

% Differential equation parameters
k = 2;
q = 2;
b = 1;

% Degree of N(u)
d = 20; 

%Bisection Method
cleft = 0.01;
cright = 2;
while abs(cright - cleft) >= 1e-5

    cmid = 0.5*(cleft + cright);
    flag = volume(cmid,d,k,q,b);
    
    if flag == 0
       cleft = cmid; 
    else
       cright = cmid;
    end
end

% Print results
fprintf('The lower bound on the minimum speed is %f found using degree %d polynomials.\n',cmid,d)

%%
function flag = volume(c,d,k,q,b)

    %SOS Parameters
    lambda = 1000;
    eps = 1e-4;

    % SDP variables
    sdpvar u v
    
    % Auxiliary function
    [V, cV] = polynomial([u v], d);

    % S Procedure
    d2 = d + q;
    [s1, c1] = polynomial([u v], d2);
    [s2, c2] = polynomial([u v], d2);
    [s3, c3] = polynomial(u, d2);

    % Derivatives
    dVdu = jacobian(V,u);
    dVdv = jacobian(V,v);

    % Replacements
    Vv0 = replace(V,v,0); %V(u,0)

    % Constraints
    cons = [];
    cons = [cons, replace(V,[u v], [0 -eps]) == 0, replace(V, [u v], [1 0]) == 0];
    cons = [cons, sos(lambda*(dVdu*v - c*dVdv*v + b*dVdv*u*v - dVdv*(1-(u^q))*u^(k+1)) - (u^k)*V - u*(1-u)*s1 + v*s2)];
    cons = [cons, sos(-Vv0 - eps*(1-u) - u*(1-u)*s3)];
    cons = [cons, sos(s1), sos(s2), sos(s3)];

    % SOS Solver
    ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
    sol = solvesos(cons,[],ops,[cV; c1; c2; c3])

    % Return whether solvesos failed or succeeded
    flag = sol.problem;
    
end






