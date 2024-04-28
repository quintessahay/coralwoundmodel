% Heat Map for Wound Size and Polyp Distance
clear all;
close all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


params=readmatrix('finalfit.txt');
params=[params; 0; 0];

% set ranges
steps = 50; %number of steps to take
wrange = [1 30]; %min and max wound size
prange = [2*((wrange(1)/pi)^(1/2)) 2*((wrange(2)/pi)^(1/2))]; %polyp distance based on wound diameter

w = wrange(1):0.1:wrange(2);
p = 0.05:0.05:(wrange(2)/pi)^(1/2);
ansmat = NaN(length(w),length(p));
u_mat = [];

m = [];

for i = 1:length(w)
    j=1;
    while p(j)<=(w(i)/pi)^(1/2)
        params(end-1) = w(i); %set wound size
        params(end) = p(j); %set polyp distance

        fun = @(t, u) coral_ODEs_rhs(t, u, params);
        tspan = 1:50000;
        Opt = odeset('Events', @(t, u) myEvent(t,u,w(i)));

        [t,u,te,ye,ie] = ode45(fun,tspan,[0,0,0,0],Opt);
        usave = interp1(t,u(:,4),tspan(1):100:tspan(end));
        u_mat = [u_mat; usave];
        ansmat(i, j) = te(1);
        m = [m;p(j) w(i) te(1)];
        if j==length(p)
            break;
        end
        j=j+1;
    end
end


for i = 1:length(m)
    if m(i,3) >= 2000
        m(i,3) = 2000;
    end
end

figure;
scatter(m(:, 1), m(:, 2), 100, m(:, 3), "filled",'s')
hold on
scatter(1.8,11.59303667,250,'kd','filled','LineWidth',2);
scatter(1.2,23.6625,250,'k','x','LineWidth',4);
xlabel('$$max_{pd}$$ (mm)')
ylabel('$$C_{\infty}$$ (mm\textsuperscript{2})')
cb = colorbar;
set(cb,'XTick',400:400:2000,'XTickLabel',{'400','800','1200','1600','$$\geq 2000$$'},'FontSize',18)
set(gca,'FontSize',30,'FontName','Arial')


function [value, isterminal, direction] = myEvent(t, u, wsize)
    value = double(abs(u(4)-wsize)>0.01*wsize);
    isterminal = 1;
    direction = 0;
end