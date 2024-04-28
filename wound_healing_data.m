list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

circ_base=[0	22.62626
23.25581	18.08081
71.70543	17.9798
166.66667	17.47475
191.86047	16.16162
238.37209	14.0404
335.27132	13.83838
360.46512	12.22222
406.97674	11.81818
536.82171	3.23232
573.64341	2.12121
594.96124	0.90909
668.60465	0.50505
691.86047	0.40404
717.05426	0.40404
738.37209	0.40404
761.62791	0.30303
835.27132	0];
circ_mid=[0	18.22581
23.4375	16.37097
72.26563	14.19355
167.96875	12.17742
191.40625	12.01613
240.23438	11.69355
335.9375	10.64516
363.28125	13.46774
410.15625	10.16129
541.01563	10.08065
578.125	10.80645
599.60938	10.48387
673.82812	10.24194
697.26563	10
722.65625	10.16129
744.14062	6.69355
767.57812	9.75806
841.79687	7.98387
863.28125	9.51613
892.57812	7.82258];
circ_top=[0	24.93976
23.4375	23.85542
72.26563	22.28916
167.96875	20.84337
191.40625	20.72289
240.23438	20.72289
335.9375	18.31325
363.28125	17.46988
410.15625	15.78313
541.01563	8.6747
578.125	6.38554
599.60938	4.21687
673.82812	1.3253
697.26563	1.20482
722.65625	0.84337
744.14062	0.48193
767.57812	0.24096
841.79687	0.24096
863.28125	0.24096
890.625	0];
circ_mean = readmatrix('circ_mean.txt');
circ_mean(:,2) = circ_mean(end,2)-circ_mean(:,2);

lin_base=[0	9.15663
23.4375	8.19277
72.26563	7.751
167.96875	7.14859
191.40625	6.62651
240.23438	6.38554
335.9375	6.10442
363.28125	5.74297
410.15625	4.97992
541.01563	4.01606
578.125	3.65462
599.60938	2.04819
673.82812	1.40562
697.26563	1.36546
722.65625	1.36546
744.14062	1.20482
767.57812	0.68273
841.79687	1.40562
863.28125	1.36546
892.57812	1.44578
914.0625	0.40161
939.45312	0.2008];
lin_mid=[0	9.59839
23.4375	9.32103
72.26563	9.28872
167.96875	8.13943
191.40625	7.86207
240.23438	6.98638
337.89062	6.76111
363.28125	6.60454
410.15625	6.33095
541.01563	4.70538
578.125	4.06878
599.60938	2.86741
673.82812	1.83515
697.26563	1.35699
722.65625	1.20043
744.14062	1.00307
767.57812	0.52491
841.79687	0];
lin_top=[0	12.77108
23.4375	11.80723
72.26563	11.20482
167.96875	10.66265
191.40625	10.48193
240.23438	10.12048
335.9375	9.6988
363.28125	9.21687
410.15625	8.79518
541.01563	3.6747
578.125	3.01205
599.60938	2.10843
673.82812	1.0241
697.26563	0.72289
722.65625	0.54217
744.14062	0.42169
767.57812	0.3012
841.79687	0.3012
863.28125	0.24096
890.625	0.18072];
lin_mean = readmatrix('lin_mean.txt');
lin_mean(:,2) = lin_mean(end,2)-lin_mean(:,2);

figure
h2=tiledlayout(1,2); % This plots the C curve with the data overlay
    nexttile
    plot(lin_base(:,1),lin_base(:,2),'*-','LineWidth',2,'MarkerSize',12)
    hold on
    plot(lin_mid(:,1),lin_mid(:,2),'s-','LineWidth',2,'MarkerSize',12)
    plot(lin_top(:,1),lin_top(:,2),'d-','LineWidth',2,'MarkerSize',12)
    plot(lin_mean(:,1),lin_mean(:,2),'.--','LineWidth',2,'MarkerSize',12)
    ylim([0 30])
    legend({'Base','Mid','Top','Mean'},'Location','northeast','Interpreter','Latex');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
    xlabel('Hours','FontSize',30,'FontName','Arial','Interpreter','Latex');
    nexttile
    plot(circ_base(:,1),circ_base(:,2),'*-','LineWidth',2,'MarkerSize',12)
    hold on
    plot(circ_mid(:,1),circ_mid(:,2),'s-','LineWidth',2,'MarkerSize',12)
    plot(circ_top(:,1),circ_top(:,2),'d-','LineWidth',2,'MarkerSize',12)
    plot(circ_mean(:,1),circ_mean(:,2),'.--','LineWidth',2,'MarkerSize',12)
    ylim([0 30])
    legend({'Base','Mid','Top','Mean'},'Location','northeast','Interpreter','Latex');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
    xlabel('Hours','FontSize',30,'FontName','Arial','Interpreter','Latex');
    ylabel(h2,'Wound Size (mm\textsuperscript{2})','FontSize',30,'FontName','Arial','Interpreter','Latex');

figure
h2=tiledlayout(1,2); % This plots the C curve with the data overlay
    nexttile
    plot(lin_base(:,1),lin_base(:,2),'*-','LineWidth',2,'MarkerSize',12)
    hold on
    plot(lin_mid(:,1),lin_mid(:,2),'s-','LineWidth',2,'MarkerSize',12)
    plot(lin_top(:,1),lin_top(:,2),'d-','LineWidth',2,'MarkerSize',12)
    %plot(lin_mean(:,1),lin_mean(:,2),'.--','LineWidth',2,'MarkerSize',12)
    ylim([0 30])
    legend({'Base','Mid','Top'},'Location','northeast','Interpreter','Latex');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
    xlabel('Hours','FontSize',30,'FontName','Arial','Interpreter','Latex');
    nexttile
    plot(circ_base(:,1),circ_base(:,2),'*-','LineWidth',2,'MarkerSize',12)
    hold on
    plot(circ_mid(:,1),circ_mid(:,2),'s-','LineWidth',2,'MarkerSize',12)
    plot(circ_top(:,1),circ_top(:,2),'d-','LineWidth',2,'MarkerSize',12)
    %plot(circ_mean(:,1),circ_mean(:,2),'.--','LineWidth',2,'MarkerSize',12)
    ylim([0 30])
    legend({'Base','Mid','Top'},'Location','northeast','Interpreter','Latex');
    set(gca,'LineWidth',2,'FontSize',30,'FontName','Arial')
    xlabel('Hours','FontSize',30,'FontName','Arial','Interpreter','Latex');
    ylabel(h2,'Wound Size (mm\textsuperscript{2})','FontSize',30,'FontName','Arial','Interpreter','Latex');
    

figure
plot(circ_base(:,1),circ_base(:,2),'*-','LineWidth',2)
hold on
plot(circ_mid(:,1),circ_mid(:,2),'*-','LineWidth',2)
plot(circ_top(:,1),circ_top(:,2),'*-','LineWidth',2)
plot(circ_mean(:,1),circ_mean(:,2),'*--','LineWidth',2)
ylim([0 25])
set(gca,'FontSize',30)
title("Circular Wound")
xlabel('Hours')
ylabel('Wound Size (mm\textsuperscript{2})')
legend({'Base','Mid','Top','Mean'},'location','northeast')


figure
plot(lin_base(:,1),lin_base(:,2),'*-','LineWidth',2)
hold on
plot(lin_mid(:,1),lin_mid(:,2),'*-','LineWidth',2)
plot(lin_top(:,1),lin_top(:,2),'*-','LineWidth',2)
plot(lin_mean(:,1),lin_mean(:,2),'*--','LineWidth',2)
ylim([0 25])
set(gca,'FontSize',30)
title("Linear Wound")
xlabel('Hours')
ylabel('Wound Size (mm\textsuperscript{2})')
legend({'Base','Mid','Top','Mean'},'location','northeast')