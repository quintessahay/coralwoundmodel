function [a0, f0] = steadystate_analytical(params)

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

a = -kaf*kaa*(kfa+muf)-mua*kaa*(kfa+muf)+kfa*kaf*kaa;
b = sc*kaa*(kfa+muf)-kaf*(musc+ka+kf)*(kfa+muf)-mua*(musc+ka+kf)*(kfa+muf)+kfa*kaf*(musc+ka+kf);
c = sc*ka*(kfa+muf)+kfa*sc*kf;

a0 = (-b-sqrt(b^2-4*a*c))/(2*a);

f0 = (sc*kf)/((musc+ka+kaa*a0+kf)*(kfa+muf))+(kaf*a0)/(kfa+muf);
end