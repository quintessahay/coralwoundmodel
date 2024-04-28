function [x,fval,exitflag,output]=cpe_main(x0)


global time1 time2 data1 data2 R1 R2 fc sol1 sol2 cinf1 cinf2 Rvec1 Rvec2 maxpd1 maxpd2

opts = optimset('fminsearch');
opts.MaxFunEvals = 100000;
opts.MaxIter = 100000;
opts.TolFun = 1e-5;
opts.TolX = 1e-5;
opts.PlotFcns = 'optimplotfval';
LB=zeros(length(x0),1);
% (1) km (2) kma (3) mum (4) sc (5) kam
% (6) kaa (7) musc (8) kaf (9) minf (10) kfa
% (11) mua (12) kf (13) muf (14) kc (15) xc 
% (16) pf (17) ka];
A = diag(ones(17,1));
b = 100*ones(17,1);
    
[x,fval,exitflag,output] = fmincon(@cpe_solver,x0,A,b,[],[],LB,[],[],opts);
end