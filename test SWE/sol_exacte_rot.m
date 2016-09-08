function [vt,rotvt] = sol_exacte_rot(x,y,z)
global teta0 teta1 u0
%% test 1

[n1,n2]=size(x);
vt=zeros(n1,n2,3);
rotvt=zeros(n1,n2,3);
[lambda, teta,rr]=cart2sph(x,y,z);

% uu=u0*fun10(teta,teta0,teta1);
% duu=u0*dfun10(teta,teta0,teta1);

pwr=10;
uu=u0*cos(teta).^pwr;
duu=-u0*pwr.*sin(teta).*cos(teta).^(pwr-1);

vv=zeros(n1,n2);

elambda(1:n1,1:n2,1) = -sin(lambda);
elambda(1:n1,1:n2,2) =  cos(lambda);
elambda(1:n1,1:n2,3)=zeros(size(x));

eteta(1:n1,1:n2,1) = -sin(teta).*cos(lambda);
eteta(1:n1,1:n2,2) = -sin(teta).*sin(lambda);
eteta(1:n1,1:n2,3) =  cos(teta);

er(1:n1,1:n2,1) = cos(teta).*cos(lambda);
er(1:n1,1:n2,2) = cos(teta).*sin(lambda);
er(1:n1,1:n2,3) = sin(teta);

for k=1:3
    vt(:,:,k)=uu.*elambda(1:n1,1:n2,k) + vv.*eteta(1:n1,1:n2,k);
end

rotvt(:,:,1)=-(duu./rr).*cos(teta).*cos(lambda);
rotvt(:,:,2)=-(duu./rr).*cos(teta).*sin(lambda);
rotvt(:,:,3)=(uu./(rr.*cos(teta)))-(duu./rr).*sin(teta);



%% test 2
% [n1,n2]=size(x);
% for i=1:n1
%     for j=1:n2
%         vt(i,j,1)=x(i,j);
%         vt(i,j,2)=y(i,j);
%         vt(i,j,3)=z(i,j);
%     end
% end
% rotvt=zeros(n1,n2,3);

% %% test 3
% [n1,n2]=size(x);
% pwr=5;
% for i=1:n1
%     for j=1:n2
%         vt(i,j,1)=(pwr+1)*x(i,j)^pwr;
%         vt(i,j,2)=(pwr+1)*y(i,j)^pwr;
%         vt(i,j,3)=(pwr+1)*z(i,j)^pwr;
%     end
% end
% rotvt=zeros(n1,n2,3);




end