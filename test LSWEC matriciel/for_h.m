function [forh] = for_h(x,y,z,t)
% forcaque sur l'équation de conservation
global test
global teta0 teta1
if test == 1
     sigma=10^-4;
    [aaa, teta,aaa]=cart2sph(x,y,z);
    forh = -sigma.*fun10(teta,teta0,teta1).*exp(-sigma*t); 

else
    forh=zeros(size(x));
end

end