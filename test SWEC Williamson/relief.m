function [hs] = relief(x,y,z)
global test

if test == 0
    hs=zeros(size(x));
    
elseif test == 1
    hs0=2000;
    R=pi/9;
    lambdac=0*pi/2;
    tetac=pi/6;
    
    [lambda, teta, ~]=cart2sph(x,y,z);
    [n1,n2]=size(x);
    hs=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            r=sqrt(min([R^2, (lambda(i,j)-lambdac)^2+(teta(i,j)-tetac)^2]));
            hs(i,j)=hs0.*(1-r/R);
        end
    end
end
end

