function [dpmm] = derlegendre(m,l,z)
% ********************************************
% Derivative of Legendre Polynomial
%
% Author :
%     - Matthieu Brachet
%
% inspired by W. H. Press and al. in 1992
%
% ********************************************;

fact=prod(1:2:2*m-1);
dpmm=0;
if (l == m)
    dpmm=-m*fact*z.*(1-z.^2)^(m/2-1);
else
    dpmmp1=-m*fact*z.*(1-z.^2)^(m/2-1);
    if (l == m+1)
        pmmp1= legendre(m,l,z);
        dpmm=(2*m-1).*pmmp1+z.*(2*m+1).*dpmmp1;
    else
        for ll=m+2:l
            p=legendre(m,ll-1,z);
            dpmm=(2*ll-1).*p+z.*(2*ll-1).*dpmm-(ll+m-1).*dpmmp1;
            dpmm=dpmm./(l-m);
        end
    end
end

end