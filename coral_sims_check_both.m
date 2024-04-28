% Check data fits using range - circular wounds
% input-matrix of params 
% output-variable transients for M A F C
clear all

input_filename='LHSmatrix.txt'; %Change for each batch of pameters

LHSmatrix=readmatrix(input_filename); % or LHSmatrix_fit_lognorm_old

if min(min(LHSmatrix))<0
    error("matrix contains negative parameters");
end

%load data
data1=readmatrix('circ_range_reduced.txt'); %max pd 1.2
data2=readmatrix('range_lin_reduced.txt'); %max pd 1.8
f=@coral_ODEs_rhs; %Differerntial equation file

tspan=0:1000; % time vector for simulation output
rows_output=length(tspan);
[runs,npars]=size(LHSmatrix);

data_flag1=zeros(runs,length(data1)); %columns indicate which data is satisfied
data_flag1(:,1)=1;
data_flag2=zeros(runs,length(data2)); %columns indicate which data is satisfied
data_flag2(:,1)=1;

tic

for param_counter=1:runs
    
    % Fit circular wound
    
    pars = LHSmatrix(param_counter,:);
    pars(18) = data1(end,3);
    pars(19) = 1.2;
    
    [a0,f0] = steadystate_analytical(pars);
    
    [B,warning_flag]=ode23s_2(@(t,y) f(t,y,pars),tspan,[0 0 0 0],[]); %using f from above. Calling built in solver ode23s (change to not stop if numerically error)
        
    t = cell2mat(B(1,1)); %t is taking out time from matrix B that ode23s_2 outputs
    y = cell2mat(B(1,2)); %has all variables from B above
    A=[t y];% [time y]
    
    if (warning_flag==1) || (min(min(A))<0)  %check if numerical errors, if so store -1's
        
        T1(:,param_counter)=-1*ones(rows_output,1);
        M1(:,param_counter)=-1*ones(rows_output,1);
        Am1(:,param_counter)=-1*ones(rows_output,1);
        F1(:,param_counter)=-1*ones(rows_output,1);
        C1(:,param_counter)=-1*ones(rows_output,1);
        
    elseif a0>1 || f0>1 || a0/f0>1 %if a0 or f0 are above threshold, store -2
      
        T1(:,param_counter)=-2*ones(rows_output,1);
        M1(:,param_counter)=-2*ones(rows_output,1);
        Am1(:,param_counter)=-2*ones(rows_output,1);
        F1(:,param_counter)=-2*ones(rows_output,1);
        C1(:,param_counter)=-2*ones(rows_output,1);
        
    else
        T1(:,param_counter)=A(:,1);
        M1(:,param_counter)=A(:,2);
        Am1(:,param_counter)=A(:,3);
        F1(:,param_counter)=A(:,4);
        C1(:,param_counter)=A(:,5);
        
        t_check1=data1(:,1)+1;
        
        data_check1=C1(t_check1,param_counter);
        
        for i=1:length(data_check1)
            if data_check1(i)>=data1(i,2) && data_check1(i)<=data1(i,3)
                data_flag1(param_counter,i)=1;
            end
        end
        
    end
    
    pars(18) = data2(end,3);
    pars(19) = 1.8;
    
    [a0,f0] = steadystate_analytical(pars);
    
    
    [B,warning_flag]=ode23s_2(@(t,y) f(t,y,pars),tspan,[0 0 0 0],[]); %using f from above. Calling built in solver ode23s (change to not stop if numerically error)
        
    t = cell2mat(B(1,1)); %t is taking out time from matrix B that ode23s_2 outputs
    y = cell2mat(B(1,2)); %has all variables from B above
    A=[t y];% [time y]
    
    if (warning_flag==1) || (min(min(A))<0)  %check if numerical errors, if so store -1's
        
        T2(:,param_counter)=-1*ones(rows_output,1);
        M2(:,param_counter)=-1*ones(rows_output,1);
        Am2(:,param_counter)=-1*ones(rows_output,1);
        F2(:,param_counter)=-1*ones(rows_output,1);
        C2(:,param_counter)=-1*ones(rows_output,1);
        
    elseif a0>1 || f0>1 || a0/f0>1 %if a0 or f0 are above threshold, store -3
        
        T2(:,param_counter)=-2*ones(rows_output,1);
        M2(:,param_counter)=-2*ones(rows_output,1);
        Am2(:,param_counter)=-2*ones(rows_output,1);
        F2(:,param_counter)=-2*ones(rows_output,1);
        C2(:,param_counter)=-2*ones(rows_output,1);
        

    else
        T2(:,param_counter)=A(:,1);
        M2(:,param_counter)=A(:,2);
        Am2(:,param_counter)=A(:,3);
        F2(:,param_counter)=A(:,4);
        C2(:,param_counter)=A(:,5);
        
        t_check2=data2(:,1)+1;
        
        data_check2=C2(t_check2,param_counter);
        
        for i=1:length(data_check2)
            if data_check2(i)>=data2(i,2) && data_check2(i)<=data2(i,3)
                data_flag2(param_counter,i)=1;
            end
        end
        
    end
    
    if mod(param_counter,500)==0
        param_counter
    end
end   
  
toc

%%
for i=1:length(data_flag1)
    sum_flags1(i,1) = sum(data_flag1(i,:));
end

for i=1:length(data_flag2)
    sum_flags2(i,1) = sum(data_flag2(i,:));
end

[val1,row1]=max(sum_flags1)
[val2,row2]=max(sum_flags2)

sum_all = [sum_flags1 sum_flags2];

%% 
% find_set = sum_all(sum_all(:,1)>8,:);
final_set_ind = sum_all(:,1)>6 & sum_all(:,2)>6;
sparse(final_set_ind)
test = sum_all(final_set_ind,:);
final_set_inds = sum_all(:,1)>=7 & sum_all(:,2)>=7;
finalset_pars = LHSmatrix(final_set_inds,:);
[nsets,~] = size(finalset_pars);
circ_mean = readmatrix('circ_data.txt');
lin_mean = readmatrix('lin_data.txt');
sse_finalpars = sum((C1(circ_mean(:,1)+1,final_set_inds)-repmat(circ_mean(:,2),1,nsets)).^2,1)+sum((C2(lin_mean(:,1)+1,final_set_inds)-repmat(lin_mean(:,2),1,nsets)).^2,1);
[val, finalset_ind]=min(sse_finalpars)
writematrix(finalset_pars(finalset_ind,:),'LHS_fit.txt')
% final_params = LHSmatrix(final_set_ind,:);
% writematrix(final_params, 'best_fits_101222.txt');
% IC_check1 = IC_ss1(final_set_ind,:);
% IC_check2 = IC_ss2(final_set_ind,:);
% IC_check1 = IC_ss1(final_set_ind,2:3);
% writematrix(IC_check1,'ICs_best_fit_101222.txt');