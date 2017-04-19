function [forv] = for_v(x,y,z,t)
global test
global gp hp radius omega
global teta0 teta1

if test == 1
    sigma=10^-5;
    [lambda, teta, ~]=cart2sph(x,y,z);
    
    uu=-sigma.*(sqrt(gp*hp)/10).*fun10(teta,teta0,teta1).*exp(-sigma*t);
    vv=(gp./radius).*dfun10(teta,teta0,teta1).*exp(-sigma*t) + 2.*omega.*sin(teta).*(sqrt(gp*hp)/10).*fun10(teta,teta0,teta1).*exp(-sigma*t);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    forv(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    forv(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    forv(:,:,3)=uu.*elambda_z + vv.*eteta_z;
    
    
else 
    [n1,n2]=size(x);
    forv=zeros(n1,n2,3);
end
end

