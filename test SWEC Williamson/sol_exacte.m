function [ht,vt] = sol_exacte(x,y,z,t)
% exact solution for Williamson test.
global test
global gp u0 h0 radius omega
global alpha

if test == 0
    % test 2 of Williamson & al.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;

elseif test == 1
    % test 5 of Williamson & al.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;
    
elseif test == 2
    % test 5 of Williamson & al. wxith  smooth mountain.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;
    
end

