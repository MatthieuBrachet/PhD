function [v]=fun4_b(x,y,z,t)
global coef
global alphad u0;
global radius
% test de Williamson
global tetac lambdac
% test de Nair et Machenhauer
global gamma rho0
% test de Nair et Jablonowski
global teta0 lambda0
% test de Salesak
global lambda1 teta1
global lambda2 teta2

if coef == 0
    [n1, n2] = size(x);
    v=zeros(n1,n2);
    aa=radius;       % dans mod53, rayon de la sphère
    h0=1000;         % hauteur BUMP
    RR=aa/3;         % parametre BUMP +- rayon du bump
    omega0=(2*pi)/(12*24*3600);    % vitesse angulaire du BUMP en RD/sec
    u0=aa*omega0;                  % module de la vitesse du BUMP a la 
                                   % surface de la terre en m/s; 
    % passage dans le repère en rotation
    Pa=zeros(3,3);Pai=zeros(3,3);
    Pa=[cos(alphad),0,-sin(alphad);0,1,0;sin(alphad),0,cos(alphad)];
    Pai=inv(Pa);
    % rotatin due au mouvement
    Rmt=zeros(3,3);
    Rmt=[cos(-omega0*t),-sin(-omega0*t),0;...
      sin(-omega0*t),cos(-omega0*t),0;
      0,0,1];
    % application sur les coordonnées
    xx=zeros(3,n1,n2);
    xx(1,1:n1,1:n2)=x;xx(2,1:n1,1:n2)=y;xx(3,1:n1,1:n2)=z;
    xxt=zeros(3,n1,n2);
    for i=1:n1,
         for j=1:n2,
        xxt(1:3,i,j)=Pa*Rmt*Pai*xx(1:3,i,j);
        end
    end
    xt=zeros(n1,n2);yt=zeros(n1,n2);zt=zeros(n1,n2);
    xt=xxt(1,1:n1,1:n2);
    yt=xxt(2,1:n1,1:n2);
    zt=xxt(3,1:n1,1:n2);
    [lambdat,tetat,~]=cart2sph(xt,yt,zt);
    % calcul de la nouvelle solution au temps t
    rdt=aa*acos(sin(tetac)*sin(tetat)+cos(tetac)*(cos(tetat).*cos(lambdat-lambdac)));
    lwkt=(rdt<RR);
    v=zeros(n1,n2);
    for i=1:n1,
        for j=1:n2,  
        v(i,j)=0.5*h0*(1+cos(pi*rdt(1,i,j)/RR))*lwkt(1,i,j);
        end
    end
    
elseif coef == 1
    
    [lambda, teta, ~]=cart2sph(x,y,z);
    [ lambda_prime, teta_prime ] = rotated_coord( lambda, teta );
    
    rho=rho0*cos(teta_prime);
    V=u0*(3*sqrt(3)/2).*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    wr=V./(radius*rho).*logi;
    v=1-tanh((rho./gamma).*sin(lambda_prime - wr*t));

elseif coef == 2 

    ws=u0/radius;
    [lambda,teta,~]=cart2sph(x,y,z);
       
%% lambda_s teta_s
     [lambda_prime, teta_prime]=rotated_coord(lambda, teta);
     lambdas_prime = lambda_prime - ws*t;
     tetas_prime=teta_prime;
     [lambdas, tetas]=unrotated_coord(lambdas_prime, tetas_prime);
     lambdatc=lambda0;
     tetatc=teta0;
    
%% solution
    nom=cos(tetas).*sin(lambdas-lambdatc);
    denom=cos(tetas).*sin(tetatc).*cos(lambdas-lambdatc)-cos(tetatc).*sin(tetas);
    lambdas_prime=atanp(nom,denom);
    tetas_prime=asin(sin(tetas).*sin(tetatc)+cos(tetas).*cos(tetatc).*cos(lambdas-lambdatc));

    rho=rho0.*cos(tetas_prime);
    V=u0*(3*sqrt(3)/2)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    wr=V./(radius.*rho).*logi;
    v=1-tanh((rho./gamma).*sin(lambdas_prime - wr*t)); 
    
elseif coef == 3
    
    [n1, n2] = size(x);
    v=zeros(n1,n2);
    h0=1000;         % hauteur BUMP
    RR=radius/3;         % parametre BUMP +- rayon du bump
    omega0=(2*pi)/(12*24*3600);    % vitesse angulaire du BUMP en RD/sec
    u0=radius*omega0;                  % module de la vitesse du BUMP a la 
                                   % surface de la terre en m/s; 
    % passage dans le repère en rotation
    Pa=zeros(3,3);Pai=zeros(3,3);
    Pa=[cos(alphad),0,-sin(alphad);0,1,0;sin(alphad),0,cos(alphad)];
    Pai=inv(Pa);
    % rotatin due au mouvement
    Rmt=zeros(3,3);
    Rmt=[cos(-omega0*t),-sin(-omega0*t),0;...
      sin(-omega0*t),cos(-omega0*t),0;
      0,0,1];
    % application sur les coordonnées
    xx=zeros(3,n1,n2);
    xx(1,1:n1,1:n2)=x;xx(2,1:n1,1:n2)=y;xx(3,1:n1,1:n2)=z;
    xxt=zeros(3,n1,n2);
    for i=1:n1,
        for j=1:n2,
            xxt(1:3,i,j)=Pa*Rmt*Pai*xx(1:3,i,j);
        end
    end
    xt=zeros(n1,n2);yt=zeros(n1,n2);zt=zeros(n1,n2);
    xt=squeeze(xxt(1,1:n1,1:n2));
    yt=squeeze(xxt(2,1:n1,1:n2));
    zt=squeeze(xxt(3,1:n1,1:n2));
    [lambdat,tetat,~]=cart2sph(xt,yt,zt);
    v=zalesak(lambdat, tetat);
    
    
    
end
end