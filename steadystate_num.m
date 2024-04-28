%This version of steadystate_num does not use ode23s_2

% This file calculates the numerical steady state for A and F
% to obtain the initial conditions for those variables
function [a0,f0]=steadystate_num(params)
    % Set the debris accumulation rate to 0 since we simulate
    % with no wound
    params(1)=0;
    
    % Define our ode file
    f=@coral_ODEs_rhs;
    
    % Define time vector for determining steady state, tolerance
    % level, and steps before to check for steady state
    tspan_ss=0:800;
    tol=0.01;
    steps_before=20;
    
    % Define two initial conditions to simulate to steady state
    % M(0)=0, A(0)=F(0)=0.5 or A(0)=F(0)=1, and C(0)=C_infinity
    y0=zeros(1,4);
    y0(2)=0;
    y0(3)=0;
    y0(4)=params(end-1);
    y1=zeros(1,4);
    y1(2)=0.5;
    y1(3)=0.5;
    y1(4)=params(end-1);
    
    % Initialize keep_ss which is a 0 or 1 variable for if the
    % simulation reached a numerical steady state
    keep_ss=0;
    
    % Simulate the odes using ode23s_2 with the first initial
    % condition and extract the time (t) and variable vectors (y)
    [t, y]=ode23s(@(t,y) f(t,y,params),tspan_ss,y0,[]);
    A=[t y];
    
    % Check if the simulation reached a numerical steady state
    % (no warning flag, no negative values, and solution is not
    % changing at the end values)
    if (min(min(A))>=-0.001) && norm(A(end-steps_before:end,2:end)-A(end,2:end))<tol
        keep_ss=1;
    end
    
    % If steady state was not reached above rerun with the 
    % second initial condition
    if keep_ss==0
        [t, y]=ode23s(@(t,y) f(t,y,params),tspan_ss,y1,[]);
        A=[t y];
        
        % Check if steady state was reached
        if (min(min(A))>=-0.001) && norm(A(end-steps_before:end,2:end)-A(end,2:end))<tol %numerical ss and no warning signs and no negative variable values
            keep_ss=1;
        end
        
    end
    
    % If numerical error occured or negative values just set
    % initial conditions of A and F to 0
    if (min(min(A))<-.001)
        a0=0;
        f0=0; 
    % If steady state was not reached just set the initial
    % conditions of A and F to 0
    elseif  keep_ss==0 
        a0=0;
        f0=0;
    % Otherwise use the end values as the intial conditions
    else
        a0=A(end,3);
        f0=A(end,4);
    end
end