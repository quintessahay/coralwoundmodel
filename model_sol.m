%--------------------------------------------------------------------------
%Solve the ODE using the given parameters (pars) and initial conditions
%(Init) at the given points (xdata)
%Set ydata to global
%In this case, x and xdata are the same so don't need to interpolate ydata
%--------------------------------------------------------------------------

function [sol,rout] = model_sol(pars,x,Init,set,ydata)

global ODE_TOL 

options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@(t,y,pars) modelBasic(t,y,pars,set),[0 x(end)],Init,options,pars);
% Evaluates solution at time-points where data are measured
sol   = deval(sol,x);

% Specify variable(s) you want sensitivity for
C = sol(4,:);
rout = C;

% Sensitivity wrt residual or scaled residual
% rout = (sol(4,:)' - ydata);
%rout = (sol(2,:)' - ydata)/mean(ydata);
