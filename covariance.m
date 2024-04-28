
%--------------------------------------------------------------------------
%Using sensitivity matrix computed in DriverBasic_sens.m, identifies 
%pairwise correlations between parameters and returns sensitive and 
%linearly independent parameter subset
%
%--------------------------------------------------------------------------
function INDMAP = covariance
%clear all; close all;

%global DIFF_INC

% DIFF_INC = 1e-4; % sqrt (ODE_TOL)

load Sens.mat %sensitivity matrix (sens), sensitivities (Rsens), and 
              %ranking (Isens)
% load sens_scale

DIFF_INC = 0.05;
Is = find(Rsens/max(Rsens)>DIFF_INC);

INDMAP = Isens(Is); %(Is)

labels = {'k_m','k_{ma}','\mu_m','s_c',...
    'k_{am}','k_{aa}','\mu_{sc}','k_{af}',...
    'm_\infty','k_{fa}','\mu_a','k_f',...
    '\mu_f','k_c','x_c','p_f','k_a'};
reduced_labels = labels(INDMAP);

S = sens(:,INDMAP);
[m,n] = size(S);

% Model Hessian
A  = S'*S;
Ai = inv(A);

disp('condition number of A = transpose(S)S and S');
disp([ cond(A) cond(S) cond(Ai)] );

% Calculate the covariance matrix
[a,b] = size(Ai);
for i = 1:a
    for j = 1:b
        r(i,j)=Ai(i,j)/sqrt(Ai(i,i)*Ai(j,j)); % covariance matrix
    end;
end;

rn = triu(r,1) % extract upper triangular part of the matrix
[i,j] = find(abs(rn)>0.95); % parameters with a value bigger than 0.95 are correlated

disp('correlated parameters');
for k = 1:length(i)
   disp([INDMAP(i(k)),INDMAP(j(k)),rn(i(k),j(k))]);
end

while ~isempty(i)
    INDMAP = INDMAP(find(INDMAP~=INDMAP(j(k))))
    pause;
    r = [];
    S = sens(:,INDMAP);

    [m,n] = size(S);
    A  = S'*S;
    size(A)
    Ai = inv(A);
    
    disp('condition number of A = transpose(S)S and S');
    disp([ cond(A) cond(S)  cond(Ai)] );

    [a,b] = size(Ai);
    
    for i = 1:a
        for j = 1:b
            r(i,j)=Ai(i,j)/sqrt(Ai(i,i)*Ai(j,j)); % covariance matrix
        end;
    end;
    
    rn = triu(r,1); % extract upper triangular part of the matrix
    [i,j] = find(abs(rn)>0.95); % parameters with a value bigger than 0.95 are correlated

    disp('correlated parameters');
    for k = 1:length(i)
       disp([INDMAP(i(k)),INDMAP(j(k)),rn(i(k),j(k))]);
    end
end
