clear all
close all

% parameter values: this file contains optimized parameters
% the only difference between the two wounds is the C_infinity
% (initial wound size - mm^2) and max_pd (maximum distance 
% between polyps
params_mat=readmatrix('final_fit.txt');
params=params_mat;
params1=params;
params2=params;


% data: these files contain real time (in hours) and closure
% data for wound healing in corals
data1=readmatrix('lin_data.txt'); %linear wound
params1(18)=data1(end,2); % cinf for linear wound
params1(19)=1.8; % max_pd for linear wound
data2=readmatrix('circ_data.txt'); %circular wound
params2(18)=data2(end,2); %cinf for circular wound
params2(19)=1.2; % max_pd for linear wound

% Initial conditions
ics = [0 0 0 0];

% Establish a time vector
tspan=(0:1:1000)'; % 0 to 1000 hours in increments of 1

% Solve the system: this solves the system of odes using the 
% built in ode23s solver in Matlab - the inputs are the ode
% file, the time vector, the initial conditions, options (blank)
% and the parameter values
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol1 = ode23s(@coral_ODEs_rhs,tspan,ics,options,params1);
t1 = sol1.x;
u1 = sol1.y;
cplot1 = u1(4,:);

sol2 = ode23s(@coral_ODEs_rhs,tspan,ics,options,params2);
t2 = sol2.x;
u2 = sol2.y;
cplot2 = u2(4,:);

solnew1=deval(sol1,data1(:,1));
R1 = solnew1(4,2:end)'-data1(2:end,2);
solnew2=deval(sol2,data2(:,1));
R2 = solnew2(4,2:end)'-data2(2:end,2);
fc_lin = R1'*R1/length(R1);
fc_circ = R2'*R2/length(R2);

%% Plot some figures:

figure
h1=tiledlayout(1,2); % This plots the C curve with the data overlay
    nexttile
    plot(t1,cplot1,'LineWidth',2);
    hold on
    plot(data1(:,1),data1(:,2),'*','LineWidth',2)
    ylabel('C (mm\textsuperscript{2})','FontSize',30,'FontName','Arial');  
    xlabel('Hours','FontSize',30,'FontName','Arial');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
    nexttile
    plot(t2,cplot2,'LineWidth',2);
    hold on
    plot(data2(:,1),data2(:,2),'*','LineWidth',2)
    xlabel('Hours','FontSize',30,'FontName','Arial');
    legend({'Model', 'Data'},'Location','southeast');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
 

resids = [R1; R2];
resid_lim = ceil(max(abs(resids)));

figure
    s1=scatter(data1(2:end,1),R1,'k','filled');
    s1.SizeData = 100;
    hold on
    s2=scatter(data2(2:end,1),R2,'k','^','LineWidth',2);
    s2.SizeData = 100;
    yline(0,'k','LineWidth',2)
    ylim([-resid_lim resid_lim])
    xlabel('Hours')
    ylabel('Residual')
    set(gca,'LineWidth',2,'FontSize',30, 'FontName','Arial')
    [hLg, icons]=legend({'','Linear','Circular'},'Location','northeast');
    icons = findobj(icons,'Type','patch');
    icons = findobj(icons,'Marker','none','-xor');
    set(icons(1),'MarkerSize',11,'LineWidth',2);

figure
h3=tiledlayout(2,2); % This plots all the model variables
    nexttile
    plot(t1,scaled_C1,'k-','LineWidth',2)
    hold on
    plot(t2,scaled_C2,'k--','LineWidth',2)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Arial')
    legend({'Linear','Circular'},'Location','southeast','FontSize',20)
    ylabel('C (normalized)');
    nexttile
    plot(t1,u1(2,:),'k-','LineWidth',2);
    hold on
    plot(t2,u2(2,:),'k--','LineWidth',2)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Arial')
    ylabel('A-units');
    nexttile
    plot(t1,u1(1,:),'k-','LineWidth',2);
    hold on
    plot(t2,u2(1,:),'k--','LineWidth',2);
    ylim([0 max(ceil(max(sol1.y(1,:))),ceil(max(sol2.y(1,:))))+1]);
    ylabel('M-units')
    xlabel('Hours')
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Arial')
    nexttile
    plot(t1,u1(3,:),'k-','LineWidth',2);
    hold on
    plot(t2,u2(3,:),'k--','LineWidth',2);
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Arial')
    ylabel('F-units');
    xlabel('Hours')
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Arial')  
    h3.TileSpacing = 'compact';
    h3.Padding = 'compact';
