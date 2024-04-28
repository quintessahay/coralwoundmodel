%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set 
%pars = log(pars) in DriverBasic_sens
%Plots ranked sensitivities 
%--------------------------------------------------------------------------

function DriverBasic_sens

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

% data1 and data2 because I have two sets of data

global ydata1 xdata1 ydata2 xdata2 

data1 = readmatrix('lin_data.txt');
data2 = readmatrix('circ_data.txt');

xdata1 = data1(2:end,1);
ydata1 = data1(2:end,2);
xdata2 = data2(2:end,1);
ydata2 = data2(2:end,2);

%senseq finds the non-weighted sensitivities
[pars,Init] = load_global;
[sens1,~,sol1] = senseq(pars,1:1000,Init,1,ydata1);
[sens2,~,sol2] = senseq(pars,1:1000,Init,2,ydata2);

sens = [sens1; sens2];
sol = [sol1(4,:)'; sol2(4,:)'];
for i=1:length(sol)
    if sol(i)==0
        sol(i)=1;
    end
end
% ranked classical sensitivities
[M,N] = size(sens);
for i=1:N
    sens(:,i)=sens(:,i).*(pars(i)./(sol));
end
for i = 1:N
  sens_norm(i)=norm(sens(:,i),2)*sqrt(1/M);
end

[Rsens,Isens] = sort(sens_norm,'descend');
display([Isens]);
par_names = {'$$k_m$$','$$k_{ma}$$','$$\mu_m$$','$$s_c$$',...
    '$$k_{am}$$','$$k_{aa}$$','$$\mu_{sc}$$','$$k_{af}$$',...
    '$$m_\infty$$','$$k_{fa}$$','$$\mu_a$$','$$k_f$$',...
    '$$\mu_f$$','$$k_c$$','$$x_c$$','$$p_f$$','$$k_a$$'};

Rsens_scaled = Rsens./max(Rsens);

%Ranked sensitivities
figure(1);clf; 
% h=semilogy(Rsens_scaled,'x');
h=plot(Rsens_scaled,'x');
hold on
yline(0.05,'--k','LineWidth',4)
set(h,'linewidth',4);
set(h,'Markersize',24);
set(gca,'Fontsize',24);
set(gca,'defaultAxesTickLabelInterpreter','latex')
grid on;
ylabel('Normalized Sensitivity Values');
xlabel('Parameters')
xticks(1:length(par_names))
xticklabels(par_names(Isens))
print -depsc2 RankedSensitivities.eps

save Sens.mat sens Rsens Isens sens_norm Rsens_scaled;