clear all
%close all
global time1 time2 data1 data2 R1 R2 fc sol1 sol2 cinf1 cinf2 Rvec1 Rvec2 maxpd1 maxpd2
Rvec=[];
datamat1=readmatrix('lin_data.txt');
datamat2=readmatrix('circ_data.txt');
time1=datamat1(:,1);
data1=datamat1(:,2);
cinf1=datamat1(end,2);
time2=datamat2(:,1);
data2=datamat2(:,2);
cinf2=datamat2(end,2);
maxpd1=1.8;
maxpd2=1.2;

params = readmatrix('LHS_bestfit.txt');
km=params(1); % rate that foreign bacteria and debris/damaged cells accumulates in the wound
kma=params(2); % rate that amoebocytes remove foreign bacteria and debris/damaged cells
mum=params(3); % decay rate of foreign bacteria/ other debris
sc=params(4); % source rate of stem cells that can become amoebocytes or fibroblasts
kam=params(5); % rate that stem cells differentiate to amoebocytes (activated by M)
kaa=params(6); % rate that stem cells differentiate to amoebocytes (activated by A)
musc=params(7); % decay rate of stem cells
kaf=params(8); % rate that amoebocytes transition to fibroblasts
minf=params(9); % value that controls inhibition from M
kfa=params(10); % rate that fibroblasts transition to amoebocytes

mua=params(11); % decay rate of amoebocytes
kf=params(12); % baseline rate that stem cells differentiate to fibroblasts
muf=params(13); % decay rate of fibroblasts
kc=params(14); % rate of epidermal coral cells maturation (activated by F)
xc=params(15); % value that controls the rate that epidermal coral cells mature (MM)
pf=params(16); % rate that fibroblasts proliferate
ka=params(17);

x0=[km; kma; mum; sc; kam; kaa; musc; kaf; minf; kfa; mua; kf; muf; kc; xc; pf; ka];
[x,fval,exitflag,output]=cpe_main(x0)

fc_lin = R1(2:end)'*R1(2:end)/length(R1(2:end));
fc_circ = R2(2:end)'*R2(2:end)/length(R2(2:end));
%%
figure
plot(sol1.x,sol1.y(4,:),time1,data1,'*','LineWidth',2)
title('Linear Wound')
subtitle(strcat(['MSE=', num2str(fc_lin)]));
legend({'Model','Data'},'Location','northwest')
xlabel('Hours')
ylabel('Amount Healed (mm^2)')
set(gca,'LineWidth',2,'FontSize',18)

figure
scatter(time1(2:end),R1,'filled')
yline(0)
title('Residual Plot')


figure
plot(sol2.x,sol2.y(4,:),time2,data2,'*','LineWidth',2)
title('Circular Wound')
subtitle(strcat(['MSE=', num2str(fc_circ)]));
legend({'Model','Data'},'Location','northwest')
xlabel('Hours')
ylabel('Amount Healed (mm^2)')
set(gca,'LineWidth',2,'FontSize',18)

figure
scatter(time2(2:end),R2,'filled')
yline(0)
title('Residual Plot')

figure
plot(sol1.x,sol1.y(2,:),'r','LineWidth',2)
hold on
plot(sol1.x,sol1.y(3,:),'b','LineWidth',2)
title('Linear Wound')
legend({'A','F'},'Location','northwest')
xlabel('Hours')
ylabel('Cells')
set(gca,'LineWidth',2,'FontSize',18)

figure
plot(sol2.x,sol2.y(2,:),'r','LineWidth',2)
hold on
plot(sol2.x,sol2.y(3,:),'b','LineWidth',2)
title('Circular Wound')
legend({'A','F'},'Location','northwest')
xlabel('Hours')
ylabel('Cells')
set(gca,'LineWidth',2,'FontSize',18)

writematrix(x,'cpe_fit.txt')
save('cpe_out.mat','fc','fc_circ','fc_lin','output','Rvec1','Rvec2','x','x0')


