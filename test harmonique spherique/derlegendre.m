function [dpmm] = derlegendre(m,l,z)
% ********************************************
% Derivative of Legendre Polynomial
%
% Author :
%     - Matthieu Brachet
%
% inspired by W. H. Press and al. in 1992
%
% ********************************************
pmm=1;
if (m>=0)
    fact=1;
    for i=1:m
        fact=fact+2;
    end
    dpmm=m.*fact.*z.*(1-z.^2).^(m/2-1);
    if (l == m)
        dpmm=dpmm;
    else
        pmm=legendre(m,m,z);
        pmmp1=dpmm;
        if (l == m+1)
            dpmm=(2*m+1)*pmm+z.*(2*m+1).*dpmm;
        else
            for ll=m+2:l
                pmm=legendre(m,ll,z);
                pll=((2*ll-1)*pmm+z.*(2*l-1).*dpmm-(ll+m-1)*pmmp1)/(ll-m);
                dpmm=pmmp1;
                pmmp1=pll;
            end
            dpmm=pll;
        end
    end
end
end