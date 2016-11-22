function [hs] = relief(x,y,z)
global test
global hs0_mount R_mount lambdac_mount tetac_mount
if test == 0
    hs=zeros(size(x));
    
elseif test == 1
    
    [lambda, teta, ~]=cart2sph(x,y,z);
    [n1,n2]=size(x);
    hs=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            r=sqrt(min([R_mount^2, (lambda(i,j)-lambdac_mount)^2+(teta(i,j)-tetac_mount)^2]));
            hs(i,j)=hs0_mount.*(1-r./R_mount);
        end
    end
    
elseif test == 2
    
    [lambda, teta, ~]=cart2sph(x,y,z);
    [n1,n2]=size(x);
    hs=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            r=sqrt(min([R_mount^2, (lambda(i,j)-lambdac_mount)^2+(teta(i,j)-tetac_mount)^2]));
            hs(i,j)=hs0_mount.*exp(-(2.8*r./R_mount).^2);
        end
    end
end

