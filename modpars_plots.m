% Heat Map for M Dynamics - km and kma
clear all;
% close all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

params=readmatrix('finalfit.txt');
params=[params; 11.59303667; 1.8];

par1 = 1;
par2 = 2;

par_names = {'$$k_m$$','$$k_{ma}$$','$$\mu_m$$','$$s_c$$',...
    '$$k_{am}$$','$$k_{aa}$$','$$\mu_{sc}$$','$$k_{af}$$',...
    '$$m_\infty$$','$$k_{fa}$$','$$\mu_a$$','$$k_f$$',...
    '$$\mu_f$$','$$k_c$$','$$x_c$$','$$p_f$$','$$k_a$$'};

% set ranges
per_adj = [1 0.5 1 1.5 1];
per_adj1 = [1 1 0.5 1 1.5];

M = [];
A = [];
F = [];
C = [];

for i = 1:length(per_adj)
    j=i;
        params_adj = params;
        params_adj(par1) = params(par1)*per_adj(i);
        params_adj(par2) = params(par2)*per_adj1(j);
        
        fun = @(t, u) coral_ODEs_rhs(t, u, params_adj);
        tspan = 1:5000;
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        [t,u] = ode23s(fun,tspan,[0 0 0 0],options);
        
        M = [M u(:,1)];
        A = [A u(:,2)];
        F = [F u(:,3)];
        C = [C u(:,4)];
end

figure
h=tiledlayout(2,2);
nexttile
plot(tspan,C(:,1),'k','LineWidth',2)
hold on
plot(tspan,C(:,2),'r--','LineWidth',2)
plot(tspan,C(:,3),'r:','LineWidth',2)
plot(tspan,C(:,4),'b--','LineWidth',2)
plot(tspan,C(:,5),'b:','LineWidth',2,'MarkerIndices',1:400:length(tspan),'MarkerSize',10)
ylabel('C-units','Interpreter','Latex')
set(gca,'FontSize',30,'FontName','Arial')
xticks([0 1000 2000])
xticklabels({'0','1','2'})
xlim([0 2500])
legend({'Baseline','Decreased $$k_{m}$$','Decreased $$k_{ma}$$','Increased $$k_{m}$$','Increased $$k_{ma}$$'},'Location','southeast','Interpreter','Latex','FontSize',20)
nexttile
plot(tspan,A(:,1),'k','LineWidth',2)
hold on
plot(tspan,A(:,2),'r--','LineWidth',2)
plot(tspan,A(:,3),'r:','LineWidth',2)
plot(tspan,A(:,4),'b--','LineWidth',2)
plot(tspan,A(:,5),'b:','LineWidth',2,'MarkerIndices',1:400:length(tspan),'MarkerSize',10)
ylabel('A-units','Interpreter','Latex')
set(gca,'FontSize',30,'FontName','Arial')
xticks([0 1000 2000])
xticklabels({'0','1','2'})
xlim([0 2500])
ylim([0 0.2])
nexttile
plot(tspan,M(:,1),'k','LineWidth',2)
hold on
plot(tspan,M(:,2),'r--','LineWidth',2)
plot(tspan,M(:,3),'r:','LineWidth',2)
plot(tspan,M(:,4),'b--','LineWidth',2)
plot(tspan,M(:,5),'b:','LineWidth',2,'MarkerIndices',1:400:length(tspan),'MarkerSize',10)
ylabel('M-units','Interpreter','Latex')
xlabel('Hours (thousands)')
set(gca,'FontSize',30,'FontName','Arial')
xticks([0 1000 2000])
xticklabels({'0','1','2'})
xlim([0 2500])
nexttile
plot(tspan,F(:,1),'k','LineWidth',2)
hold on
plot(tspan,F(:,2),'r--','LineWidth',2)
plot(tspan,F(:,3),'r:','LineWidth',2)
plot(tspan,F(:,4),'b--','LineWidth',2)
plot(tspan,F(:,5),'b:','LineWidth',2,'MarkerIndices',1:400:length(tspan),'MarkerSize',10)
ylabel('F-units','Interpreter','Latex')
xlabel('Hours (thousands)')
set(gca,'FontSize',30,'FontName','Arial')
xticks([0 1000 2000])
xticklabels({'0','1','2'})
xlim([0 2500])
h.TileSpacing = 'compact';
h.Padding = 'compact';

%% Increase/decrease both
% Heat Map for AF Dynamics - s_c
clear all;
% close all;

params=readmatrix('finalfit_041024.txt');
params=[params; 11.59303667; 1.8];

par1 = 4;

par_names = {'$$k_m$$','$$k_{ma}$$','$$\mu_m$$','$$s_c$$',...
    '$$k_{am}$$','$$k_{aa}$$','$$\mu_{sc}$$','$$k_{af}$$',...
    '$$m_\infty$$','$$k_{fa}$$','$$\mu_a$$','$$k_f$$',...
    '$$\mu_f$$','$$k_c$$','$$x_c$$','$$p_f$$','$$k_a$$'};

% set ranges
per_adj = [1 0.9 1.1];

M = [];
A = [];
F = [];
C = [];

for i = 1:length(per_adj)
        params_adj = params;
        params_adj(par1) = params(par1)*per_adj(i);
        fun = @(t, u) coral_ODEs_rhs(t, u, params_adj);
        tspan = 1:5000;
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        [t,u] = ode23s(fun,tspan,[0 0 0 0],options);
        
        M = [M u(:,1)];
        A = [A u(:,2)];
        F = [F u(:,3)];
        C = [C u(:,4)];
end

figure;
h=tiledlayout(2,2);
nexttile
plot(tspan,C(:,1),'k','LineWidth',3)
hold on
plot(tspan,C(:,2),'k--','LineWidth',3)
plot(tspan,C(:,3),'k:','LineWidth',3)
ylabel('C-units')
set(gca,'FontSize',30,'FontName','Arial')
xlim([0 2000])
xticks([0 1000 2000])
xticklabels({'0','1','2'})
legend({'Baseline','10\% Decrease','10\% Increase'},'Location','southeast','FontSize',20)
nexttile
plot(tspan,A(:,1),'k','LineWidth',3)
hold on
plot(tspan,A(:,2),'k--','LineWidth',3)
plot(tspan,A(:,3),'k:','LineWidth',3)
ylabel('A-units')
set(gca,'FontSize',30,'FontName','Arial')
xlim([0 2000])
ylim([0 0.2])
xticks([0 1000 2000])
xticklabels({'0','1','2'})
nexttile
plot(tspan,M(:,1),'k','LineWidth',3)
hold on
plot(tspan,M(:,2),'k--','LineWidth',3)
plot(tspan,M(:,3),'k:','LineWidth',3)
ylabel('M-units')
xlabel('Hours (thousands)')
set(gca,'FontSize',30,'FontName','Arial')
xlim([0 2000])
xticks([0 1000 2000])
xticklabels({'0','1','2'})
nexttile
plot(tspan,F(:,1),'k','LineWidth',3)
hold on
plot(tspan,F(:,2),'k--','LineWidth',3)
plot(tspan,F(:,3),'k:','LineWidth',3)
ylabel('F-units')
xlabel('Hours (thousands)')
set(gca,'FontSize',30,'FontName','Arial')
xlim([0 2000])
xticks([0 1000 2000])
xticklabels({'0','1','2'})
h.TileSpacing = 'compact';
h.Padding = 'compact';
