function [pmm] = legendre(m,l,x)
% ********************************************
% Legendre Polynomial
%
% Author :
%     - Matthieu Brachet
%
% inspired by W. H. Press and al. in 1992
%
% ********************************************
pmm=1;
if (m>=0)
    somx2=sqrt((1-x).*(1+x));
    fact=1;
    for i=1:m
        pmm=fact*somx2*pmm;
        fact=fact+2;
    end
    if (l == m)
        pmm=pmm;
    else
        pmmp1=x*(2*m+1)*pmm;
        if (l == m+1)
            pmm=pmmp1;
        else
            for ll=m+2:l
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm=pmmp1;
                pmmp1=pll;
            end
            pmm=pll;
        end
    end
end
end