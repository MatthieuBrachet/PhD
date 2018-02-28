function [ht] = sol_exacte(x,y,z,t)
global test
global radius

if test == 0
    %% test 3 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    xx=x./radius; yy=y./radius; zz=z./radius;
    ht=(xx+yy+zz)./sqrt(3);
elseif test == 1
    %% test 1 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    [lambda, teta, ~]=cart2sph(x,y,z);
    ht=sin(lambda-pi).*(teta>-pi/12).*(teta<pi/12);
elseif test == 2
    %% test 4 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    [lambda, teta, ~]=cart2sph(x,y,z);
    ht=fun1(x).*cos(lambda).*cos(teta);
end  
end

