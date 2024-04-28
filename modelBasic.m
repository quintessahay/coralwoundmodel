%--------------------------------------------------------------------------
%Define ODEs, using parameters from load_global.m 
%--------------------------------------------------------------------------

function ydot = modelBasic(t,y,pars,set)

M  = y(1);
A  = y(2);
F  = y(3);
C  = y(4);

km=pars(1); % rate that foreign bacteria and debris/damaged cells accumulates in the wound
kma=pars(2); % rate that amoebocytes remove foreign bacteria and debris/damaged cells
mum=pars(3); % decay rate of foreign bacteria/ other debris
sc=pars(4); % source rate of stem cells that can become amoebocytes or fibroblasts
kam=pars(5); % rate that stem cells differentiate to amoebocytes (activated by M)
kaa=pars(6); % rate that stem cells differentiate to amoebocytes (activated by A)
musc=pars(7); % INSENSITIVE decay rate of stem cells
kaf=pars(8); % rate that amoebocytes transition to fibroblasts
minf=pars(9); % value that controls inhibition from M
kfa=pars(10); % rate that fibroblasts transition to amoebocytes

mua=pars(11); % decay rate of amoebocytes
kf=pars(12); % INSENSITIVE rate that stem cells differentiate to fibroblasts
muf=pars(13); % decay rate of fibroblasts
kc=pars(14); % rate of epidermal coral cells maturation (activated by F)
xc=pars(15); % value that controls the rate that epidermal coral cells mature (MM)
pf=pars(16); % fibroblast proliferation from other fibroblasts and lack on contact inhibition
ka=pars(17); % INSENSITIVE baseline rate at which cells differentiate to amoebocytes
if set==1
    cinf=11.59303667; % initial wound size (mm^2)
    maxpd=1.8; % maximum distance between polyps in the wound (mm)
else
    cinf=23.511015; % initial wound size (mm^2)
    maxpd=1.2; % maximum distance between polyps in the wound (mm)
end

dy1=km*(1-C/cinf)-kma*M*A-mum*M; %debris

dy2=(sc*(1/(1+(maxpd*(1-C/cinf))))*(ka+kam*M+kaa*A))/(musc+ka+kam*M+kaa*A+(kf+pf*F*(1-C/cinf)))-kaf*A*(1/(1+(M/minf)))+kfa*F-mua*A; %amoebocytes

dy3=(sc*(1/(1+(maxpd*(1-C/cinf))))*(kf+pf*F*(1-C/cinf)))/(musc+(kf+pf*F*(1-C/cinf))+ka+kam*M+kaa*A)+kaf*A*(1/(1+(M/minf)))-kfa*F-muf*F; %fibroblasts

dy4=(kc*F)/(xc+F)*(1-C/cinf); %coral epidermal cells

ydot = [dy1; dy2; dy3; dy4];