function [v,dv]=fun6(x,y,z)
% % n=max(size(x,1),size(x,2));
% % v=zeros(n,1);
% % v1=zeros(n,1);v2=zeros(n,1);v3=zeros(n,1);
% %
% n1=size(x,1);
% n2=size(x,2);
% v=zeros(n1,n2,3);
% dv=zeros(n1,n2);
% % % --1--
% % % FORMULA MBA/JF BASD ON A SCALAR PRODUCT.
% % % THIS GIVES A NULL DIV VECTOR.
% % % VECTOR FIELD F(X)=N(X) \CROSS PHI
% % phi=zeros(n1,n2,3);
% % %
% % %phi(1:n1,1:n2,1)=0;phi(1:n1,1:n2,2)=0;phi(1:n1,1:n2,3)=0;
% % phi(1:n1,1:n2,1)=4;phi(1:n1,1:n2,2)=1;phi(1:n1,1:n2,3)=3;
% % %
% % xnorm=zeros(n1,n2,3);
% % nnorm=sqrt(x.*x+y.*y+z.*z);
% % xnorm(1:n1,1:n2,1)=x./nnorm;
% % xnorm(1:n1,1:n2,2)=y./nnorm;
% % xnorm(1:n1,1:n2,3)=z./nnorm;
% % v(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,3)...
% %     - xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,2);
% % v(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,1)...
% %     - xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,3);
% % v(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,2)...
% %     - xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,1);
% % dv(1:n1,1:n2)=zeros(n1,n2);
% % % --2--
% % % VECTOR FIELD F(X)=C(X) (N(X) \CROSS PHI)
% % % DIV_S(F)=GRAD_S C(X) \DOT (N(X) \CROSS PHI)
% phi=zeros(n1,n2,3);
% cv=zeros(n1,n2);
% gcv=zeros(n1,n2,3);
% %---------------------------------------------
% phi(1:n1,1:n2,1)=1;phi(1:n1,1:n2,2)=2;phi(1:n1,1:n2,3)=3;
% % phi(1:n1,1:n2,1)=1;phi(1:n1,1:n2,2)=2;phi(1:n1,1:n2,3)=3;
% %----------------------
% a=1;
%  kw1=1;kw2=2;kw3=3;
% %  kw1=2;kw2=6;kw3=10;
% %kw1=1;kw2=1;kw3=1;
% cv=a*(sin(kw1*pi*x)+sin(kw2*pi*y)+sin(kw3*pi*z));
% gcv(1:n1,1:n2,1)=a*kw1*pi*cos(kw1*pi*x);
% gcv(1:n1,1:n2,2)=a*kw2*pi*cos(kw2*pi*y);
% gcv(1:n1,1:n2,3)=a*kw3*pi*cos(kw3*pi*z);
% %
% xnorm=zeros(n1,n2,3);
% nnorm=sqrt(x.*x+y.*y+z.*z);
% xnorm(1:n1,1:n2,1)=x./nnorm;
% xnorm(1:n1,1:n2,2)=y./nnorm;
% xnorm(1:n1,1:n2,3)=z./nnorm;
% %
% gscv(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)...
%                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
% gscv(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)...
%                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
% gscv(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)...
%                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);
% %               
% v1=zeros(n1,n2,3);
% v1(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,3)...
%     - xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,2);
% v1(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,1)...
%     - xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,3);
% v1(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,2)...
%     - xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,1);
% %
% v=zeros(n1,n2,3); % CHAMP DE VECTEURS (OUT)
% v(1:n1,1:n2,1)=cv(1:n1,1:n2).*v1(1:n1,1:n2,1);
% v(1:n1,1:n2,2)=cv(1:n1,1:n2).*v1(1:n1,1:n2,2);
% v(1:n1,1:n2,3)=cv(1:n1,1:n2).*v1(1:n1,1:n2,3);
% %
% dv(1:n1,1:n2)=zeros(n1,n2); % DIVERGENCE CHAMP DE VECTEURS (OUT)
% dv(1:n1,1:n2)=dot(gscv(1:n1,1:n2,1:3),v1(1:n1,1:n2,1:3),3);
% % % --2 bis--
% % % VECTOR FIELD F(X)=C(X) (N(X) \CROSS PHI)
% % % DIV_S(F)=GRAD_S C(X) \DOT (N(X) \CROSS PHI)
% % phi=zeros(n1,n2,3);
% % cv=zeros(n1,n2);
% % gcv=zeros(n1,n2,3);
% % %---------------------------------------------
% % phi(1:n1,1:n2,1)=1;phi(1:n1,1:n2,2)=1;phi(1:n1,1:n2,3)=1;
% % %----------------------
% % a=1;
% % kw1=1;kw2=1;kw3=1;
% % cv=a*(exp(kw1*x)+exp(kw2*y)+exp(kw3*z));
% % gcv(1:n1,1:n2,1)=a*kw1*exp(kw1*x);
% % gcv(1:n1,1:n2,2)=a*kw2*exp(kw2*y);
% % gcv(1:n1,1:n2,3)=a*kw3*exp(kw3*z);
% % %
% % xnorm=zeros(n1,n2,3);
% % nnorm=sqrt(x.*x+y.*y+z.*z);
% % xnorm(1:n1,1:n2,1)=x./nnorm;
% % xnorm(1:n1,1:n2,2)=y./nnorm;
% % xnorm(1:n1,1:n2,3)=z./nnorm;
% % %
% % gscv(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)...
% %                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
% % gscv(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)...
% %                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
% % gscv(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)...
% %                    -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);
% % %               
% % v1=zeros(n1,n2,3);
% % v1(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,3)...
% %     - xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,2);
% % v1(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,1)...
% %     - xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,3);
% % v1(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,2)...
% %     - xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,1);
% % %
% % v=zeros(n1,n2,3); % CHAMP DE VECTEURS (OUT)
% % v(1:n1,1:n2,1)=cv(1:n1,1:n2).*v1(1:n1,1:n2,1);
% % v(1:n1,1:n2,2)=cv(1:n1,1:n2).*v1(1:n1,1:n2,2);
% % v(1:n1,1:n2,3)=cv(1:n1,1:n2).*v1(1:n1,1:n2,3);
% % %
% % dv(1:n1,1:n2)=zeros(n1,n2); % DIVERGENCE CHAMP DE VECTEURS (OUT)
% % dv(1:n1,1:n2)=dot(gscv(1:n1,1:n2,1:3),v1(1:n1,1:n2,1:3),3);
% % % --3--
% % VECTOR FIELD F(X)=(N(X) \CROSS GRAD H(x))
% % DIV_S(F)=0, CF PAPER JCP DE M. BEN-ARTZI ET AL.(GEOMETRU COMPATIBLE
% % TANGENTIAL VECTOR FIELDS)
% % hh=zeros(n1,n2);
% % ghh=zeros(n1,n2,3);
% % 
% %    a=1;
% %    kw1=2;kw2=6;kw3=10;
% %    hh=a*(sin(kw1*pi*x)+sin(kw2*pi*y)+sin(kw3*pi*z));
% %    ghh(1:n1,1:n2,1)=a*kw1*pi*cos(kw1*pi*x);
% %    ghh(1:n1,1:n2,2)=a*kw2*pi*cos(kw2*pi*y);
% %    ghh(1:n1,1:n2,3)=a*kw3*pi*cos(kw3*pi*z);
% %  
% % xnorm=zeros(n1,n2,3);
% % nnorm=sqrt(x.*x+y.*y+z.*z);
% % xnorm(1:n1,1:n2,1)=x./nnorm;
% % xnorm(1:n1,1:n2,2)=y./nnorm;
% % xnorm(1:n1,1:n2,3)=z./nnorm;
% % 
% % v=zeros(n1,n2,3);% CHAMP DE VECTEURS (OUT)
% % v(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*ghh(1:n1,1:n2,3)...
% %     - xnorm(1:n1,1:n2,3).*ghh(1:n1,1:n2,2);
% % v(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*ghh(1:n1,1:n2,1)...
% %     - xnorm(1:n1,1:n2,1).*ghh(1:n1,1:n2,3);
% % v(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*ghh(1:n1,1:n2,2)...
% %     - xnorm(1:n1,1:n2,2).*ghh(1:n1,1:n2,1);
% % 
% %  dv(1:n1,1:n2)=zeros(n1,n2); % DIVERGENCE CHAMP DE VECTEURS (OUT)
% % % --4--
% % VECTOR FIELD F(X)=(N(X) \CROSS GRAD H(x))
% % DIV_S(F)=0, CF PAPER JCP DE M. BEN-ARTZI ET AL.(GEOMETRU COMPATIBLE
% % TANGENTIAL VECTOR FIELDS)
% % hh=zeros(n1,n2);
% % ghh=zeros(n1,n2,3);
% % % 
% %    a=1;
% %    kw1=1;kw2=1;kw3=1;
% %    hh=a*(exp(kw1*x)+exp(kw2*y)+exp(kw3*z));
% %    ghh(1:n1,1:n2,1)=a*kw1*exp(kw1*x);
% %    ghh(1:n1,1:n2,2)=a*kw2*exp(kw2*y);
% %    ghh(1:n1,1:n2,3)=a*kw3*exp(kw3*z);
% %  %
% % xnorm=zeros(n1,n2,3);
% % nnorm=sqrt(x.*x+y.*y+z.*z);
% % xnorm(1:n1,1:n2,1)=x./nnorm;
% % xnorm(1:n1,1:n2,2)=y./nnorm;
% % xnorm(1:n1,1:n2,3)=z./nnorm;
% % %
% % v=zeros(n1,n2,3);% CHAMP DE VECTEURS (OUT)
% % v(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*ghh(1:n1,1:n2,3)...
% %     - xnorm(1:n1,1:n2,3).*ghh(1:n1,1:n2,2);
% % v(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*ghh(1:n1,1:n2,1)...
% %     - xnorm(1:n1,1:n2,1).*ghh(1:n1,1:n2,3);
% % v(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*ghh(1:n1,1:n2,2)...
% %     - xnorm(1:n1,1:n2,2).*ghh(1:n1,1:n2,1);
% % %
% % dv(1:n1,1:n2)=zeros(n1,n2); % DIVERGENCE CHAMP DE VECTEURS (OUT)
% % =========================================================================
% % =========================================================================
% % =========================================================================
% % --3--
% %  v=exp(5*x)+y.^4+z.^11;
% %   v1=5*exp(5*x);
% %   v2=4*y.^3;
% %   v3=11*z.^10;
% % --4--
% %      kw1=1;kw2=2;kw3=4;
% %       v=sin(2*kw1*x*pi)+sin(2*kw2*y*pi)+sin(2*kw3*z*pi);
% %       v1=2*kw1*pi*cos(2*kw1*x*pi);
% %       v2=2*kw2*pi*cos(2*kw2*y*pi);
% %       v3=2*kw3*pi*cos(2*kw3*z*pi);
% % --5--
%   % a=6.37122d+06;
% %   a=1;
% %   kw1=1;kw2=1;kw3=1;
% %    v=a*(exp(kw1*x)+exp(kw2*y)+exp(kw3*z));
% %     v1=a*kw1*exp(kw1*x);
% %     v2=a*kw2*exp(kw2*y);
% %     v3=a*kw3*exp(kw3*z);
% % --6--
% %     kw1=1;kw2=8;kw3=15;
% %     v=exp(cos(kw1*x)).*exp(cos(kw2*y)).*exp(cos(kw3*z));
% %     v1=-kw1*sin(kw1*x).*exp(cos(kw1*x)).*exp(cos(kw2*y)).*exp(cos(kw3*z));
% %     v2=-kw2*sin(kw2*y).*exp(cos(kw2*y)).*exp(cos(kw1*x)).*exp(cos(kw3*z));
% %     v3=-kw3*sin(kw3*z).*exp(cos(kw3*z)).*exp(cos(kw1*x)).*exp(cos(kw2*y));
% % --7-- 
% % v=z;
% % v1=0*x;
% % v2=0*y;
% % v3=ones(n1,n2);
% % --8--
% % nhs=14; mhs=9; % number of the normalized spherical harmonics, pbar^m_n (Rokhlin-Tygert)
% % %teta=zeros(n1,n2);lambda=zeros(n1,n2);radius=zeros(n1,n2);
% % [lambda,teta,radius]=cart2sph(x,y,z);
% % %size(teta)
% % ic=complex(0,1);
% % %teta1=zeros(n1*n2);
% %  nwk1=zeros(n1,n2);
% % %for j=1:n2,
% % vc=zeros(n1,n2);
% % for i=1:n1,
% %     for j=1:n2,
% %     nwk=legendre(nhs,sin(teta(i,j)),'norm');
% %     %size(nwk)
% %     nwk1(i,j)=nwk(mhs+1);
% %     vc(i,j)=nwk1(i,j)*exp(ic*mhs*lambda(i,j));
% %     end
% % end
% %nwk=legendre(nhs,-sin(teta),'norm'); % REGARDER EXPLICATIONS DE CETTE ROUTINE DANS MATLAB;
% %disp('START');
% %size(nwk)
% %size(teta)
% %[s1,s2,s3]=size(nwk)
% %nwk1=zeros(size(teta));
% % nwk1=squeeze(nwk(mhs+1,:,:)); % suppression de la taille "1" en debut de tableau;
% %nwk1=squeeze(nwk(mhs+1,:)); % suppression de la taille "1" en debut de tableau;
% %size(nwk1)
% %size(nwk(mhs+1,:,:))
% %nwk1=zeros(s2,s3);
% % for i=1:s2,
% %     for j=1:s3,
% %         nwk1(i,j)=nwk(mhs+1,i,j);
% %     end
% % end
% %nwk1=nwk(mhs+1,1:s2,1:s3);
% %nwk1=squeeze(nwk1);
% %size(nwk1)
% %size(nwk(mhs+1,:,:))
% %size(lambda)
% % for i=1:n1,
% %  for j=1:n2,
% %   v(i,j)=nwk(mhs+1,i,j)*exp(ic*mhs*lambda(i,j));
% %  end
% % end
% %v=nwk1.*exp(ic*mhs*lambda);
% % for i=1:n1,
% %  for j=1:n2,
% %   v(i,j)=nwk(mhs+1,i,j)*exp(ic*mhs*lambda(i,j));% V COMPLEXE
% %  end
% % end
% % v=real(vc);
% % GRADIENT NON AFFECTES POUR L'INSTANT.
%     


%% vitesse 1:3
    [n1,n2]=size(x);
    v=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=cos( teta );
    vv=zeros(n1,n2);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    v(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    v(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    v(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    dv=zeros(size(x));
