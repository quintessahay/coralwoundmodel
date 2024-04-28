% Check collinearity combinations 4
clear all
close all

%   4     1    11    16 -- From correlation.m

col_ind = [];
load Sens.mat sens Isens Rsens
DIFF_INC = 0.05;
Is = find(Rsens/max(Rsens)>DIFF_INC);
INDMAP = Isens(Is); %(Is)
S = sens(:,INDMAP);
[~,m]=size(S);
param_vec = INDMAP;

labels = {'k_m','k_{ma}','\mu_m','s_c',...
    'k_{am}','k_{aa}','\mu_{sc}','k_{af}',...
    'm_\infty','k_{fa}','\mu_a','k_f',...
    '\mu_f','k_c','x_c','p_f','k_a'};
reduced_labels = cell2mat(labels(INDMAP));

set_size = 4;

for j=set_size
choose_mat = nchoosek(param_vec,j);
[N,~] = size(choose_mat);
CI = zeros(N,1);

for i=1:N
    sorted = sort(choose_mat(i,:));
    sub_mat = sens(:,sorted)'*sens(:,sorted);
    min_eigs = eig(sub_mat);
    min_eig=min(min_eigs);
    CI(i,1) = 1/sqrt(min_eig);
    RMS(i,1) = sum(Rsens(:,sorted));
    col_ind = [col_ind; j CI(i,1) RMS(i,1)];
end

[min_CI(j), ind(j)] = min(CI);

choose_mat(ind(j),:)
end

FIM = S'*S;
% Cond = 1.7069e+11
eigFIM=eig(FIM);
eigFIM(:,2)=eigFIM(:,1)./(max(eigFIM(:,1)));

figure
plot(log(flip(eigFIM(:,1))),'*-','LineWidth',2)
hold on
yline(0,'--','LineWidth',2)
set(gca,'FontSize',18,'FontName','Arial')
xlabel('Rank')
ylabel('Eigenvalue (Log Scale)')


C_mat = inv(FIM);
Cor_mat = zeros(size(C_mat));
for i=1:size(C_mat)
    for j=1:size(C_mat)
        Cor_mat(i,j) = C_mat(i,j)/sqrt(C_mat(i,i)*C_mat(j,j));
    end
end


low_CI = [choose_mat(col_ind(:,2)<20,:) col_ind(col_ind(:,2)<20,2) RMS(col_ind(:,2)<20,1)];
low_CI = sortrows(low_CI,set_size+1);

for i=1:length(low_CI)
    flag = 0;
    for j=2:4
        if Cor_mat(find(INDMAP==low_CI_sc(i,1)),find(INDMAP==low_CI_sc(i,j)))>0.95
            flag = 1;
        end
    end
    for j=3:4
        if Cor_mat(find(INDMAP==low_CI_sc(i,2)),find(INDMAP==low_CI_sc(i,j)))>0.95
            flag = 1;
        end
    end
        
uncorr_set_flag(i) = flag;
end

uncorr_sets = low_CI(uncorr_set_flag==0,:);

for i=1:length(uncorr_sets)
    param_uncorr_labels(i,:) = [labels(uncorr_sets(i,1:set_size)) uncorr_sets(i,set_size+1) uncorr_sets(i,set_size+2)];
end
