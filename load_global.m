%--------------------------------------------------------------------------
%Set parameters, initial conditions, and upper and lower bounds of ODE
%model
%--------------------------------------------------------------------------

function [pars,Init] = load_global
global ODE_TOL DIFF_INC

ODE_TOL  = 1e-8;
DIFF_INC = sqrt(ODE_TOL);

% Parameters file
pars = readmatrix('cpe_fit.txt');

% Initial conditions
[a0,f0] = steadystate_analytical(pars);

Init = [0 0 0 0];


