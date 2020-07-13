%--------------------------------------------------------------------------
% Application of the volume method for upper bounds on the minimum wave 
% speed to the two-component reaction-diffusion equation:
%
%    a_t = a_xx - ab^m,
%    b_t = Db_xx + ab^m,    D > 0, m - positive integer
%
% This code is associated to the paper "Minimum wave speeds in monostable 
% reaction-diffusion equations: sharp bounds by polynomial optimization" by 
% Jason J. Bramburger and David Goluskin (2020). This script is used to 
% create the data in Figure 3 and Table 4. 
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
m = 2;
D = 0.5;

% Bounding method parameters
lambda = 0.5;
d = 10;

%Bisection Method
cleft = 0;
cright = 1;
while abs(cright - cleft) >= 1e-5

    cmid = 0.5*(cleft + cright);
    flag = volume(cmid,d,D,m,lambda);

    if flag == 0
       cright = cmid; 
    else
       cleft = cmid;
    end
    
end

fprintf('An upper bound on the minimum speed is %f found using degree %d polynomials.\n',cmid,d)

%--------------------------------------------------------------------------
% Theoretical upper and lower bounds coming from: 
% "Sharp estimates on minimum traveling wave spped of reaction diffusion
% systems modelling autocatalysis' by X. Chen and Y. Qi (2007)

if D < 1
    c_upper = 4*D/sqrt(1 + 4*D);
    fprintf('The upper bound coming from Chen and Qi (2017) is %f.\n',c_upper)

    % Theoretical Lower Bounds (same reference as upper bounds)
    if m == 2
        c_lower = D/sqrt(2);
        fprintf('The lower bound coming from Chen and Qi (2017) is %f.\n',c_lower)
    end

elseif (D > 1) && (m == 2)
    c_upper = sqrt(D/(1 + 1/D));
    c_lower = sqrt(D/2);
    fprintf('The upper bound coming from Chen and Qi (2017) is %f.\n',c_upper)
    fprintf('The lower bound coming from Chen and Qi (2017) is %f.\n',c_lower)
end
%--------------------------------------------------------------------------

%%
function flag = volume(c,d,D,m,lambda)

    % Variables
    sdpvar u v w

    % Epsilon value
    eps = 1e-4;

    % Auxiliary function
    [V, cV] = polynomial([u v w], [d d d]);
    
    % S procedure polynomials
    d2 = d;
    [s1, c1] = polynomial([u v w], d2);
    [s2, c2] = polynomial([u v w], d2);
    [s3, c3] = polynomial([u v w], d2);
    [s4, c4] = polynomial([u v w], d2);

    % Derivatives
    dVdu = jacobian(V,u);
    dVdv = jacobian(V,v);
    dVdw = jacobian(V,w);

    if D < 1
        
        % Replacements
        Vv1 = replace(V,v,1); %V(u,1,w)

        % S procedure polynomials
        [s5, c5] = polynomial([u w], d2);
        [s6, c6] = polynomial([u w], d2);

        % Contraints
        cons = [];
        cons = [cons, replace(V, [u v w], [1 1 0]) == 0];
        cons = [cons, sos(-lambda*(dVdu*(v + w - u)*D + dVdv*w + dVdw*(-w + (D/(c^2))*u*((1-v)^m))) - V - u*(1-u)*s1 - v*(1-v)*s2 - w*s3 - (v-u)*s4)];
        cons = [cons, replace(V, [u v w], [0 0 0]) <= -eps];
        cons = [cons, sos(Vv1 - eps*w - u*(1-u)*s5 - w*s6)];
        cons = [cons, sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)];

        %SOS Solver
        ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
        sol = solvesos(cons,[],ops,[cV; c1; c2; c3; c4; c5; c6]);

    elseif D > 1
        
        % Replacements
        Vu1 = replace(V,u,1); %V(1,v,w)

        % S procedure polynomials
        [s5, c5] = polynomial([u w], d2);
        [s6, c6] = polynomial([u w], d2);

        % Contraints
        cons = [];
        cons = [cons, replace(V, [u v w], [1 1 0]) == 0];
        cons = [cons, sos(-lambda*(dVdu*(v + w - u)*D + dVdv*w + dVdw*(-w + (D/(c^2))*u*((1-v)^m))) - V - u*(1-u)*s1 - v*(1-v)*s2 - w*s3 - (u-v)*s4)];
        cons = [cons, replace(V, [u v w], [0 0 0]) <= -eps];
        cons = [cons, sos(Vu1 - eps*w - v*(1-v)*s5 - w*s6)];
        cons = [cons, sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)];

        %SOS Solver
        ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
        sol = solvesos(cons,[],ops,[cV; c1; c2; c3; c4; c5; c6]);
    end

    % Return whether solvesos failed or succeeded
    flag = sol.problem;

end







