function [u] = funfu( teta )
global test u0

%% data for this function
if test == 0
    u=cos(teta);
elseif test == 1
    gamma=pi/18;
    phi0=pi/4;
    u=u0.*cos(pi/(2*gamma).*(teta-phi0)).^2;
elseif test == 2
    phi0=pi/7;
    phi1=-phi0;
    e=exp(-4/(phi1-phi0).^2);
    denom=(teta-phi0).*(teta-phi1);
    u=u0./e.*exp(1./denom).*(teta <= phi0).*(phi1 <= teta);
end

end