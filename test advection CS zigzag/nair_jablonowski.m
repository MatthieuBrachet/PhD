function [vitx, vity, vitz] = nair_jablonowski(x,y,z,t)

global u0 alphad radius
global lambda0 teta0
global rho0

[lambda, teta, ~]=cart2sph(x,y,z);
ws=u0/radius;

%% position du centre du vortex

[lambda0_prime, teta0_prime] = rotated_coord(lambda0,teta0);
lambdac_prime=lambda0_prime+ws*t;
tetac_prime=teta0_prime;
[lambdatc, tetatc]=unrotated_coord(lambdac_prime,tetac_prime);

%% generation de wr

tetaprime=asin(sin(teta).*sin(tetatc)+cos(teta).*cos(tetatc).*cos(lambda-lambdatc));
rho=rho0.*cos(tetaprime);
V=u0*(3*sqrt(3)/2).*sech(rho).^2.*tanh(rho);
logi=(rho~=0);
wr=V./(radius.*rho).*logi;

%% calcul de la vitesse

us=u0.*(cos(teta).*cos(alphad)+sin(teta).*cos(lambda).*sin(alphad));
ur=radius.*wr.*(sin(tetatc).*cos(teta)-cos(tetatc).*cos(lambda-lambdatc).*sin(teta));
u =us+ur;

vs=-u0.*sin(lambda).*sin(alphad);
vr=radius.*wr.*(cos(tetatc).*sin(lambda-lambdatc));
v =vs+vr;

%% changement de coordonnées

elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = 0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
    
vitx=u.*elambda_x + v.*eteta_x;
vity=u.*elambda_y + v.*eteta_y;
vitz=u.*elambda_z + v.*eteta_z;

end