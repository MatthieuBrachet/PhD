function [ u ] = u_test(phi)
global test
global gamma u0 phi0
global phi1 umax 
if test == 0
    u = u0.*cos(pi./2.*(phi-phi0)./gamma).^2;

elseif test == 1
    for i=1:size(phi,1)
        for j=1:size(phi,2)
            pp=phi(i,j);
            if (pp >= phi0)
                u(i,j)=0;
            elseif (pp > phi0) && (pp<phi1)
                denom=(pp-phi0).*(pp-phi1);
                en=exp(-4./(phi1-phi0).^2);
                u(i,j)=umax./en.*exp(1./denom);
            else
                u(i,j)=0;
                
            end
        end
    end
end