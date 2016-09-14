function [vt,rotvt] = sol_exacte_rot(x,y,z)
global teta0 teta1 u0

[n1,n2]=size(x);
vt=zeros(n1,n2,3);
rotvt=zeros(n1,n2,3);
[lambda, teta,rr]=cart2sph(x,y,z);

uu=u0*fun10(teta,teta0,teta1);
duu=u0*dfun10(teta,teta0,teta1);

% pwr=10;
% uu=u0*cos(teta).^pwr;
% duu=-u0*pwr.*sin(teta).*cos(teta).^(pwr-1);


elambda(1:n1,1:n2,1) = -sin(lambda);
elambda(1:n1,1:n2,2) =  cos(lambda);
elambda(1:n1,1:n2,3)=zeros(size(x));

eteta(1:n1,1:n2,1) = -sin(teta).*cos(lambda);
eteta(1:n1,1:n2,2) = -sin(teta).*sin(lambda);
eteta(1:n1,1:n2,3) =  cos(teta);

for k=1:3
    glambda(:,:,k)=(1./(rr.*cos(teta))).*elambda(:,:,k);
    gteta(:,:,k)=(1./rr).*eteta(:,:,k);
end

for k=1:3
    vt(:,:,k)=uu.*elambda(1:n1,1:n2,k);
end

% du/dlambda
dulambda(:,:,1)=uu.*(-cos(lambda));
dulambda(:,:,2)=uu.*(-sin(lambda));
dulambda(:,:,3)=zeros(n1,n2,1);

% du/dteta
for k=1:3
    duteta(:,:,k)=duu.*elambda(:,:,k);
end
% calcul du rot
for i=1:n2
    for j=1:n2
        rotvt(i,j,1:3)=cross(glambda(i,j,1:3),dulambda(i,j,1:3))+cross(gteta(i,j,1:3),duteta(i,j,1:3));
    end
end



end