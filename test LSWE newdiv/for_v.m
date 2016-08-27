function [forv] = for_v(x,y,z,t)
global test
global gp hp radius omega

if test == 3
    sigma=10^-4;
    cste=sqrt(gp*hp)/10;

    [lambda, teta, ~]=cart2sph(x,y,z);
    %teta=teta+pi;
    % u_lambda
    uu = -sigma*cste*cos(teta).*exp(-sigma*t)+(gp./radius).*sin(lambda).*exp(-sigma*t);
    % u_teta
    vv = 2*omega*cste.*sin(teta).*cos(teta).*exp(-sigma*t)-(gp/radius).*sin(teta).*sin(lambda).*exp(-sigma*t);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    forv(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    forv(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    forv(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    
elseif test == 4
    sigma=10^-4;
    cste=sqrt(gp*hp)/10;
    pwa=1;
    pwb=5;

    [lambda, teta, ~]=cart2sph(x,y,z);
    dhdl=-3*pwa.*cos(3*lambda).^(pwa-1).*sin(3*lambda).*cos(teta).^pwb.*sin(sigma*t);
    dhdt=-pwb.*sin(teta).*cos(3*lambda).^pwa.*cos(teta).^(pwb-1).*sin(sigma*t);
    % u_lambda
    uu = sigma*cste*cos(teta).*cos(sigma*t)+(gp./radius).*dhdl./cos(teta);
    % u_teta
    vv = 2*omega*cste.*sin(teta).*cos(teta).*sin(sigma*t)+(gp/radius).*dhdt;

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

