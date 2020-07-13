%--------------------------------------------------------------------------
% Application of the surface method for upper bounds on the minimum wave 
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
cleft = 0;
cright = 2;
while abs(cright - cleft) >= 1e-6

    cmid = 0.5*(cleft + cright);
    flag = surface(cmid,d,k,q,b);
    
    if flag == 0
       cright = cmid; 
    else
       cleft = cmid;
    end
end

% Print results
fprintf('The upper bound on the minimum speed is %f found using degree %d polynomials.\n',cmid,d)

%%
function flag = surface(c,d,k,q,b)
 
    % SDP variables
    sdpvar u
    z = sdpvar(2,1);

    % Auxiliary function
    [N, cN] = polynomial(u,d,1);
    intN = int(-c*N + b*u*N,u,0,u);

    % S Procedure
    d2 = d + q;
    [s11, s11c] = polynomial(u,d2);
    [s12, s12c] = polynomial(u,d2);
    [s22, s22c] = polynomial(u,d2);
    S = [s11, s12; s12, s22]; %Symmetric 2x2 matrix of auxillary functions

    % Function defintions and integral terms
    r = 1/(q + k + 2)*u^(q+k+2) - 1/(k+2)*u^(k+2);
    Q = [(intN+r), N; N, 2];

    %Constraints
    cons = [sos(z'*(Q - u*(1-u)*S)*z), sos(z'*S*z)];

    %SOS Solver
    ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
    sol = solvesos(cons,[],ops,[cN; s11c; s12c; s22c]);

    % Return whether solvesos failed or succeeded
    flag = sol.problem;

end




