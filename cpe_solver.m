function [fc] = cpe_solver(g)

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

global time1 time2 data1 data2 R1 R2 fc sol1 sol2 cinf1 cinf2 Rvec1 Rvec2 maxpd1 maxpd2

% Weights for objective function
w1 = ones(length(data1),1);
w2 = ones(length(data2),1);

pars = g

[a0,f0]=steadystate_analytical(pars);
if a0<=0 || f0<=0
    fc=1e10;
elseif a0>f0
    fc=1e10;
else
    a0=0;
    f0=0;
    m0=0;
    c0=0;

    sol1=ode23s(@(t,y,pars) cpeode(t,y,pars),time1,[m0,a0,f0,c0],options,pars);
    solnew1=deval(sol1,time1);

    R1 = solnew1(4,2:end)'-data1(2:end);
    Rvec1 = [Rvec1; R1];

    sol2=ode23s(@(t,y,pars) cpeode1(t,y,pars),time2,[m0,a0,f0,c0],options,pars);    
    solnew2=deval(sol2,time2);

    R2 = solnew2(4,2:end)'-data2(2:end);
    Rvec2 = [Rvec2; R2];

    J = [R1; R2];

    fc = J'*J/(length(J));
end
%------

function dy = cpeode(t,y,pars)
    km=pars(1); kma=pars(2); mum=pars(3); sc=pars(4); kam=pars(5); kaa=pars(6);
    musc=pars(7); kaf=pars(8); minf=pars(9); kfa=pars(10); mua=pars(11);
    kf=pars(12); muf=pars(13); kc=pars(14); xc=pars(15); pf=pars(16); ka=pars(17);
    
    M=y(1); A=y(2); F=y(3); C=y(4);

    dy=zeros(4,1);

    dy(1)=km*(1-C/cinf1)-kma*M*A-mum*M;
    
    dy(2)=(sc*(1/(1+(maxpd1*(1-C/cinf1))))*(ka+kam*M+kaa*A))/(musc+ka+kam*M+kaa*A+(kf+pf*F*(1-C/cinf1)))-kaf*A*(1/(1+(M/minf)))+kfa*F-mua*A;
    
    dy(3)=(sc*(1/(1+(maxpd1*(1-C/cinf1))))*(kf+pf*F*(1-C/cinf1)))/(musc+(kf+pf*F*(1-C/cinf1))+ka+kam*M+kaa*A)+kaf*A*(1/(1+(M/minf)))-kfa*F-muf*F;
    
    dy(4)=(kc*F)/(xc+F)*(1-C/cinf1);
       
end

function dy = cpeode1(t,y,pars)
    km=pars(1); kma=pars(2); mum=pars(3); sc=pars(4); kam=pars(5); kaa=pars(6);
    musc=pars(7); kaf=pars(8); minf=pars(9); kfa=pars(10); mua=pars(11);
    kf=pars(12); muf=pars(13); kc=pars(14); xc=pars(15); pf=pars(16); ka=pars(17);
    
    M=y(1); A=y(2); F=y(3); C=y(4);

    dy=zeros(4,1);

    dy(1)=km*(1-C/cinf2)-kma*M*A-mum*M; % debris in wound
    
    dy(2)=(sc*(1/(1+(maxpd2*(1-C/cinf2))))*(ka+kam*M+kaa*A))/(musc+ka+kam*M+kaa*A+(kf+pf*F*(1-C/cinf2)))-kaf*A*(1/(1+(M/minf)))+kfa*F-mua*A; % amoebocytes
    
    dy(3)=(sc*(1/(1+(maxpd2*(1-C/cinf2))))*(kf+pf*F*(1-C/cinf2)))/(musc+(kf+pf*F*(1-C/cinf2))+ka+kam*M+kaa*A)+kaf*A*(1/(1+(M/minf)))-kfa*F-muf*F; % fibroblasts
    
    dy(4)=(kc*F)/(xc+F)*(1-C/cinf2); % healed coral cells 
       
end
end


