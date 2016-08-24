function [forv] = for_v(x,y,z,t)
global test
global gp radius omega

% forcaque sur l'équation de moment
    sigma=1/5000;
    
if test == 3

    [lambda, teta, ~]=cart2sph(x,y,z);
    
    % u_lambda
    uu=gp/radius.*cos(lambda).*exp(-sigma*t);
    
    % u_teta
    vv=-gp/radius.*sin(teta).*sin(lambda).*exp(-sigma*t)+2*omega*sin(teta).*cos(teta).*exp(-sigma*t);

    uu=-sigma.*cos(teta).*exp(-sigma*t)+uu;

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    forv(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    forv(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    forv(:,:,3)=uu.*elambda_z + vv.*eteta_z;

else 
    forv=zero(size(x),3);
end
end

