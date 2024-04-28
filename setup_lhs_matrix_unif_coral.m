% set up LHS matrix from scratch
% log uniform distributions
% Fall 2021 - Coral Model

npars=17;
range=[.01 100]; % all parameters will have the same range
params_range=[range(1)*ones(1,npars); range(2)*ones(1,npars)];

threshold=1; % sample on log scale
threshold_lin=1e10; % sample on linear scale

runs=100000;
distrib='unif';

LHSmatrix = LHS_Call(params_range(1,:), 0, params_range(2,:), 0, runs, distrib, threshold);

writematrix(LHSmatrix, 'LHSmatrix.txt')