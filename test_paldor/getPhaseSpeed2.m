function C=getPhaseSpeed2(waveFlag)
% parameters
omega = 7.29212e-5;          % Earth's angular frequency (rad./sec)
g     = 10.0;                % Earth's gravitational acceleration (m./sec.^2)
a     = 6371220.0;           % Earth's mean radius (m)
H0    = 5.0e3;               % Layer's mean depth (m)
nw     = 5;                   % chosen mode number (5)
kw     = 10;                  % chosen wave number (10)
sigma = 0.5+(0.25+kw.^2).^0.5;  % see Eq.10 in text
En     = g.*H0./a.^2.*(nw+sigma).^2;
delta0 = 3.*kw.^2.*En;
delta4 = -54.*g.*H0.*omega.*kw.^4./a.^2;

Cj = zeros(1,3);

for j=1:3
    deltaj = (0.5.*(delta4+sqrt(delta4.^2-4.*delta0.^3))).^(1./3).*exp(2.*pi.*1i.*j./3);
    Cj(j)  = -1./(3.*kw.^2).*deltaj.*(1+delta0./deltaj.^2);
end

Cj = real(Cj);

switch waveFlag
    case 'Rossby'
        C = -min(abs(Cj));
    case 'WIG'      
        C = min(Cj);
    case 'EIG'
        C = max(Cj);
end
end   