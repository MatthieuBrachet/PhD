function [ht,vt] = sol_exacte(x,y,z,time)
global waveFlag
lev=1;
[lambda,teta,~]=cart2sph(x,y,z);
[u,v,h,Ps,T]=paldor_init(teta,lambda,lev,time,waveFlag);

elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));

eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);

vt(:,:,1)=u.*elambda_x + v.*eteta_x;
vt(:,:,2)=u.*elambda_y + v.*eteta_y;
vt(:,:,3)=u.*elambda_z + v.*eteta_z;

ht=h;
end

