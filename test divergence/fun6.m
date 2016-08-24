function [v,dv]=fun6(x,y,z)
global test
% n=max(size(x,1),size(x,2));
% v=zeros(n,1);
% v1=zeros(n,1);v2=zeros(n,1);v3=zeros(n,1);
%
n1=size(x,1);
n2=size(x,2);
v=zeros(n1,n2,3);
dv=zeros(n1,n2);


if test ==1
    phi=zeros(n1,n2,3);
    cv=zeros(n1,n2);
    gcv=zeros(n1,n2,3);
    %---------------------------------------------
    phi(1:n1,1:n2,1)=1;phi(1:n1,1:n2,2)=2;phi(1:n1,1:n2,3)=3;
    % phi(1:n1,1:n2,1)=1;phi(1:n1,1:n2,2)=2;phi(1:n1,1:n2,3)=3;
    %----------------------
    a=1;
     kw1=1;kw2=2;kw3=3;
    %  kw1=2;kw2=6;kw3=10;
    %kw1=1;kw2=1;kw3=1;
    cv=a*(sin(kw1*pi*x)+sin(kw2*pi*y)+sin(kw3*pi*z));
    gcv(1:n1,1:n2,1)=a*kw1*pi*cos(kw1*pi*x);
    gcv(1:n1,1:n2,2)=a*kw2*pi*cos(kw2*pi*y);
    gcv(1:n1,1:n2,3)=a*kw3*pi*cos(kw3*pi*z);
    %
    xnorm=zeros(n1,n2,3);
    nnorm=sqrt(x.*x+y.*y+z.*z);
    xnorm(1:n1,1:n2,1)=x./nnorm;
    xnorm(1:n1,1:n2,2)=y./nnorm;
    xnorm(1:n1,1:n2,3)=z./nnorm;
    %
    gscv(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)...
                       -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
    gscv(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)...
                       -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
    gscv(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)...
                       -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);
    %               
    v1=zeros(n1,n2,3);
    v1(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,3)...
        - xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,2);
    v1(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,1)...
        - xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,3);
    v1(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,2)...
        - xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,1);
    %
    v=zeros(n1,n2,3); % CHAMP DE VECTEURS (OUT)
    v(1:n1,1:n2,1)=cv(1:n1,1:n2).*v1(1:n1,1:n2,1);
    v(1:n1,1:n2,2)=cv(1:n1,1:n2).*v1(1:n1,1:n2,2);
    v(1:n1,1:n2,3)=cv(1:n1,1:n2).*v1(1:n1,1:n2,3);
    %
    dv(1:n1,1:n2)=zeros(n1,n2); % DIVERGENCE CHAMP DE VECTEURS (OUT)
    dv(1:n1,1:n2)=dot(gscv(1:n1,1:n2,1:3),v1(1:n1,1:n2,1:3),3);

elseif test==2

%% vitesse 1:3
    [n1,n2]=size(x);
    v=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=80*cos( teta ).^5;
    vv=zeros(n1,n2);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    v(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    v(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    v(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    dv=zeros(size(x));

elseif test == 3
    
   [lambda, teta,r]=cart2sph(x,y,z);
    uu=80*cos(teta).^20;
    duu=-1600*sin(teta).*cos(teta).^19;
    
    v(:,:,1) = -sin(lambda).*uu;
    v(:,:,2) = zeros(size(x));
    v(:,:,3) = cos(lambda).*uu;
    
    dv_lambda=2.*cos(lambda).*sin(lambda).*uu;
    uteta=cos(lambda).*uu.*(sin(lambda).*sin(teta)+cos(teta));
    duteta=duu.*cos(lambda).*(sin(lambda).*sin(teta)+cos(teta))+cos(lambda).*uu.*(sin(lambda).*cos(teta)-sin(teta));
    dv_teta=-sin(teta).*uteta+cos(teta).*duteta;
    dv=(dv_lambda+dv_teta)./(r.*cos(teta));
    
elseif test == 4
    radius=sqrt(x.^2+y.^2+z.^2);
    lambda=atan2(-z,x);
    teta=asin(y./radius);
    
    uu=80*cos(teta).^20;
    vv=zeros(size(x));
    
    elambda_x=-sin(lambda);
    elambda_y=zeros(size(x));
    elambda_z=-cos(lambda);
    
    eteta_x=-sin(teta).*cos(lambda);
    eteta_y=cos(teta);
    eteta_z=sin(teta).*sin(lambda);
    
    v(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    v(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    v(:,:,3)=uu.*elambda_z + vv.*eteta_z;
    
    dv=zeros(size(x));
end
