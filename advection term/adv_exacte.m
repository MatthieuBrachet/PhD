function[adve]=adv_exacte(x,y,z)
global u0 teta0 teta1 radius
[lambda, teta, rr]=cart2sph(x,y,z);
uu=-u0^2/radius.*fun10(teta,teta0,teta1).*dfun10(teta,teta0,teta1).*tan(teta);
vv=zeros(size(x));

elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=zeros(size(x));
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);

adve(:,:,1)=uu.*elambda_x + vv.*eteta_x;
adve(:,:,2)=uu.*elambda_y + vv.*eteta_y;
adve(:,:,3)=uu.*elambda_z + vv.*eteta_z;