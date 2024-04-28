
function dy = coral_ODEs_rhs(t,y,params)

% Rename the parameters to their actual names
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
kf=params(12); % rate that stem cells differentiate to fibroblasts
muf=params(13); % decay rate of fibroblasts
kc=params(14); % rate of epidermal coral cells maturation (activated by F)
xc=params(15); % value that controls the rate that epidermal coral cells mature (MM)
pf=params(16); % fibroblast proliferation from other fibroblasts and lack on contact inhibition
ka=params(17); % baseline rate at which cells differentiate to amoebocytes
cinf=params(18); % initial wound size (mm^2)
maxpd=params(19); % maximum distance between polyps in the wound (mm)

% Establish our model variables
M=y(1); A=y(2); F=y(3); C=y(4);

% Initilize the dy vector for the model variables
dy=zeros(4,1);

% Define the odes
dy(1)=km*(1-C/cinf)-kma*M*A-mum*M; %debris

dy(2)=(sc*(1/(1+(maxpd*(1-C/cinf))))*(ka+kam*M+kaa*A))/(musc+ka+kam*M+kaa*A+(kf+pf*F*(1-C/cinf)))-kaf*A*(1/(1+(M/minf)))+kfa*F-mua*A; %amoebocytes

dy(3)=(sc*(1/(1+(maxpd*(1-C/cinf))))*(kf+pf*F*(1-C/cinf)))/(musc+(kf+pf*F*(1-C/cinf))+ka+kam*M+kaa*A)+kaf*A*(1/(1+(M/minf)))-kfa*F-muf*F; %fibroblasts

dy(4)=(kc*F)/(xc+F)*(1-C/cinf); %coral epithelial cells

end


