function [vt,rotvt] = sol_exacte_rot(x,y,z)
global teta0 teta1
% %% test 1
% u0=80;
% 
% [n1,n2]=size(x);
% vt=zeros(n1,n2,3);
% rotvt=zeros(n1,n2,3);
% [lambda, teta,rr]=cart2sph(x,y,z);
% 
% % uu=u0*fun10(teta,teta0,teta1);
% % duu=u0*dfun10(teta,teta0,teta1);
% 
% pwr=4;
% uu=u0*cos(teta).^pwr;
% duu=-u0*pwr.*sin(teta).*cos(teta).^(pwr-1);
% 
% vv=zeros(n1,n2);
% 
% elambda(1:n1,1:n2,1) = -sin(lambda);
% elambda(1:n1,1:n2,2) =  cos(lambda);
% elambda(1:n1,1:n2,3)=zeros(size(x));
% eteta(1:n1,1:n2,1) = -sin(teta).*cos(lambda);
% eteta(1:n1,1:n2,2) = -sin(teta).*sin(lambda);
% eteta(1:n1,1:n2,3) =  cos(teta);
% er(1:n1,1:n2,1) = cos(teta).*cos(lambda);
% er(1:n1,1:n2,2) = cos(teta).*sin(lambda);
% er(1:n1,1:n2,3) = sin(teta);
% 
% for k=1:3
%     vt(:,:,k)=uu.*elambda(1:n1,1:n2,k) + vv.*eteta(1:n1,1:n2,k);
%     rotvt(:,:,k)=(-1./rr).*(tan(teta).*uu+duu).*er(1:n1,1:n2,k)+(1./rr).*uu.*eteta(1:n1,1:n2,k);
% end
% 
%
% %% test 2
% [n1,n2]=size(x);
% for i=1:n1
%     for j=1:n2
%         %r=sqrt(x(i,j)^2+y(i,j)^2+z(i,j)^2);
%         vt(i,j,1)=x(i,j);
%         vt(i,j,2)=y(i,j);
%         vt(i,j,3)=z(i,j);
%     end
% end

%% test 3
[n1,n2]=size(x);
for i=1:n1
    for j=1:n2
        %r=sqrt(x(i,j)^2+y(i,j)^2+z(i,j)^2);
        vt(i,j,1)=4*x(i,j)^3;
        vt(i,j,2)=4*y(i,j)^3;
        vt(i,j,3)=4*z(i,j)^3;
    end
end


rotvt=zeros(n1,n2,3);




end