% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
global scheme
global opt_ftr
% -------------------------------------------------------------------------
global n nn;
global mm na nb;
global radius;
global xi eta dxi deta xx yy delta deltab dga;
global alfa beta;
global alfacr betacr;
global alfa1;
global alfag betag;
global pg kg;
global weights
% -------------------------------------------------------------------------
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
% -------------------------------------------------------------------------
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gr_I gr_II gr_III gr_IV gr_V gr_VI
% -------------------------------------------------------------------------
global ftr
% -------------------------------------------------------------------------
global pmat kmat
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp betas gamas
global m1 m2
global umat1 lmat1
% -------------------------------------------------------------------------
global hs0_mount R_mount lambdac_mount tetac_mount
global omega gp

nn=n+2;
mm=((nn-1)/2)+1;
na=4*(nn-1);
nb=na;

%% physical data
radius=6.37122d+06;
omega=7.292d-05;
gp=9.80616;

%% data for the mountain
hs0_mount=2000;
R_mount=pi/9;
lambdac_mount=-pi/2;
tetac_mount=pi/6;

%% --- Points on the CS
xi=linspace(-pi/4, pi/4, nn); 
dxi=(pi/2)/(nn-1);
eta=linspace(-pi/4, pi/4, nn); 
deta=(pi/2)/(nn-1);
xx=zeros(nn,nn);
yy=zeros(nn,nn);
for i=1:nn,
    xwk=tan(xi(i));
  for j=1:nn,
    xx(i,j)=xwk;
  end
end
for j=1:nn,
    xwk=tan(eta(j));
  for i=1:nn,
    yy(i,j)=xwk;
  end
end
delta=zeros(nn,nn);
for i=1:nn,
  for j=1:nn,
    delta(i,j)=1+xx(i,j)^2+yy(i,j)^2;
  end
end
deltab=zeros(n,n);
deltab=sqrt(delta); 
dga=zeros(nn,nn);
dga=(radius^2)*((1+xx.^2).*(1+yy.^2))./(delta.*deltab); 
alfa=zeros(nn,nn);
beta=zeros(nn,nn);
for j=1:nn,
    jbar=j-mm;
    etabar=jbar*deta;
    xwk=sqrt(1+tan(etabar)^2);
    for i=1:nn,
        alfa(i,j)=atan(tan(xi(i))/xwk);
    end
end
for i=1:nn,
    ibar=i-mm;
    xibar=ibar*dxi;
    xwk=sqrt(1+tan(xibar)^2);
    for j=1:nn,
        beta(i,j)=atan(tan(eta(j))/xwk);
    end
end
% -------------------------------------------------------------------------
alfacr=zeros(nn,nn);
betacr=zeros(nn,nn);
for j=1:nn,
    for i=1:nn,
        betacr(i,j)=atan(sin(-xi(i))*tan(eta(j)));
    end
end
for i=1:nn
    for j=1:nn,
        alfacr(i,j)=acos((cos(betacr(i,j))/cos(eta(j)))*sin(-xi(i)));
    end
end

alfa1=zeros(nn,nn); 
for i=1:nn,
    for j=1:nn,
        alfa1(i,j)=betacr(j,i);
    end
end

% -------------------------------------------------------------------------
alfag=zeros(4*(nn-1),nn); 
for j=1:nn,
  alfag(1:nn-1,j)=alfa(1:nn-1,j);
  alfag(nn:2*(nn-1),j)=alfacr(1:nn-1,j);
  alfag(2*(nn-1)+1:3*(nn-1),j)=alfa(1:nn-1,j)+pi;
  alfag(3*(nn-1)+1:4*(nn-1),j)=alfacr(1:nn-1,j)+pi;
end
betag=alfag'; 
% - FACE I -
for i=1:nn,  
    for j=1:nn,
        jbar=j-mm;
        etabar=jbar*deta;
        al1=alfa(i,j);
        x_fI(i,j)=cos(al1)*cos(etabar);
        y_fI(i,j)=sin(al1);
        z_fI(i,j)=cos(al1)*sin(etabar);
    end
end
x_fI=radius*x_fI;
y_fI=radius*y_fI;
z_fI=radius*z_fI;
% - FACE II -
for i=1:nn, 
    for j=1:nn,
        jbar=j-mm;
        etabar=jbar*deta;
        al2=alfa(i,j);
        x_fII(i,j)=-sin(al2);
        y_fII(i,j)=cos(al2)*cos(etabar);
        z_fII(i,j)=cos(al2)*sin(etabar);
    end
end
x_fII=radius*x_fII;
y_fII=radius*y_fII;
z_fII=radius*z_fII;
% - FACE III -
for i=1:nn,  
    for j=1:nn,
        jbar=j-mm;
        etabar=jbar*deta;
        al3=alfa(i,j);
        x_fIII(i,j)=-cos(al3)*cos(etabar);
        y_fIII(i,j)=-sin(al3);
        z_fIII(i,j)=cos(al3)*sin(etabar);
    end
end
x_fIII=radius*x_fIII;
y_fIII=radius*y_fIII;
z_fIII=radius*z_fIII;
% - FACE IV -
for i=1:nn,  
    for j=1:nn,
        jbar=j-mm;
        etabar=jbar*deta;
        al4=alfa(i,j);
        x_fIV(i,j)=sin(al4);
        y_fIV(i,j)=-cos(al4)*cos(etabar);
        z_fIV(i,j)=cos(al4)*sin(etabar);
    end
end
x_fIV=radius*x_fIV;
y_fIV=radius*y_fIV;
z_fIV=radius*z_fIV;
% - FACE V -
for i=1:nn,  
    for j=1:nn,
        jbar=j-mm;etabar=jbar*deta;
        al5=alfa(i,j);
        x_fV(i,j)=-cos(al5)*sin(etabar);
        y_fV(i,j)=sin(al5);
        z_fV(i,j)=cos(al5)*cos(etabar);
    end
end
x_fV=radius*x_fV;
y_fV=radius*y_fV;
z_fV=radius*z_fV;
% - FACE VI -
for i=1:nn, 
    for j=1:nn,
        jbar=j-mm;etabar=jbar*deta;
        al6=alfa(i,j);
        x_fVI(i,j)=cos(al6)*sin(etabar);
        y_fVI(i,j)=sin(al6);
        z_fVI(i,j)=-cos(al6)*cos(etabar);
    end
end
x_fVI=radius*x_fVI;
y_fVI=radius*y_fVI;
z_fVI=radius*z_fVI;
%% --- vector on the CS
% - PANEL I -
gxi_I=zeros(nn,nn,3);geta_I=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=x_fI(i,j)*(1+xx(i,j)^2);
      gxi_I(i,j,1)= -xx(i,j);
      gxi_I(i,j,2)= 1;
      gxi_I(i,j,3)= 0;
      %
      gxi_I(i,j,1:3)=gxi_I(i,j,1:3)/xwk;
    end
end

for i=1:nn,
    for j=1:nn,
      xwk=x_fI(i,j)*(1+yy(i,j)^2);
      geta_I(i,j,1)= -yy(i,j);
      geta_I(i,j,2)= 0;
      geta_I(i,j,3)= 1;
      %
      geta_I(i,j,1:3)=geta_I(i,j,1:3)/xwk;     
    end
end
% - PANEL II -
 gxi_II=zeros(nn,nn,3);geta_II=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=y_fII(i,j)*(1+xx(i,j)^2);
      gxi_II(i,j,1)= -1;
      gxi_II(i,j,2)= -xx(i,j);
      gxi_II(i,j,3)= 0;
      %
      gxi_II(i,j,1:3)=gxi_II(i,j,1:3)/xwk;
    end
end
%
for i=1:nn,
    for j=1:nn,
      xwk=y_fII(i,j)*(1+yy(i,j)^2);
      geta_II(i,j,1)= 0;
      geta_II(i,j,2)= -yy(i,j);
      geta_II(i,j,3)= 1;
      %
      geta_II(i,j,1:3)=geta_II(i,j,1:3)/xwk;     
    end
end
% - PANEL III -
 gxi_III=zeros(nn,nn,3);geta_III=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=x_fIII(i,j)*(1+xx(i,j)^2);
      gxi_III(i,j,1)= -xx(i,j);
      gxi_III(i,j,2)= 1;
      gxi_III(i,j,3)= 0;
      %
      gxi_III(i,j,1:3)=gxi_III(i,j,1:3)/xwk;
    end
end
%
for i=1:nn,
    for j=1:nn,
      xwk=x_fIII(i,j)*(1+yy(i,j)^2);
      geta_III(i,j,1)= -yy(i,j);
      geta_III(i,j,2)= 0;
      geta_III(i,j,3)= -1;
      %
      geta_III(i,j,1:3)=geta_III(i,j,1:3)/xwk;     
    end
end
% - PANEL IV -
 gxi_IV=zeros(nn,nn,3);geta_IV=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=y_fIV(i,j)*(1+xx(i,j)^2);
      gxi_IV(i,j,1)= -1;
      gxi_IV(i,j,2)= -xx(i,j);
      gxi_IV(i,j,3)= 0;
      %
      gxi_IV(i,j,1:3)=gxi_IV(i,j,1:3)/xwk;
    end
end
%
for i=1:nn,
    for j=1:nn,
      xwk=y_fIV(i,j)*(1+yy(i,j)^2);
      geta_IV(i,j,1)= 0;
      geta_IV(i,j,2)= -yy(i,j);
      geta_IV(i,j,3)= -1;
      %
      geta_IV(i,j,1:3)=geta_IV(i,j,1:3)/xwk;     
    end
end
% - PANEL V -
 gxi_V=zeros(nn,nn,3);geta_V=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=z_fV(i,j)*(1+xx(i,j)^2);
      gxi_V(i,j,1)= 0;
      gxi_V(i,j,2)= 1;
      gxi_V(i,j,3)= -xx(i,j);
      %
      gxi_V(i,j,1:3)=gxi_V(i,j,1:3)/xwk;
    end
end
%
for i=1:nn,
    for j=1:nn,
      xwk=z_fV(i,j)*(1+yy(i,j)^2);
      geta_V(i,j,1)= -1;
      geta_V(i,j,2)= 0;
      geta_V(i,j,3)= -yy(i,j);
      %
      geta_V(i,j,1:3)=geta_V(i,j,1:3)/xwk;     
    end
end
% - PANEL VI -
 gxi_VI=zeros(nn,nn,3);geta_VI=zeros(nn,nn,3); 
for i=1:nn,
    for j=1:nn,
      xwk=z_fVI(i,j)*(1+xx(i,j)^2);
      gxi_VI(i,j,1)= 0;
      gxi_VI(i,j,2)= -1;
      gxi_VI(i,j,3)= -xx(i,j);
      %
      gxi_VI(i,j,1:3)=gxi_VI(i,j,1:3)/xwk;
    end
end
%
for i=1:nn,
    for j=1:nn,
      xwk=z_fVI(i,j)*(1+yy(i,j)^2);
      geta_VI(i,j,1)= -1;
      geta_VI(i,j,2)= 0;
      geta_VI(i,j,3)= -yy(i,j);
      %
      geta_VI(i,j,1:3)=geta_VI(i,j,1:3)/xwk;     
    end
end

for i=1:nn
    for j=1:nn
        gr_I(i,j,1)=x_fI(i,j)/radius;
        gr_I(i,j,2)=y_fI(i,j)/radius;
        gr_I(i,j,3)=z_fI(i,j)/radius;
        %
        gr_II(i,j,1)=x_fII(i,j)/radius;
        gr_II(i,j,2)=y_fII(i,j)/radius;
        gr_II(i,j,3)=z_fII(i,j)/radius;
        %
        gr_III(i,j,1)=x_fIII(i,j)/radius;
        gr_III(i,j,2)=y_fIII(i,j)/radius;
        gr_III(i,j,3)=z_fIII(i,j)/radius;
        %
        gr_IV(i,j,1)=x_fIV(i,j)/radius;
        gr_IV(i,j,2)=y_fIV(i,j)/radius;
        gr_IV(i,j,3)=z_fIV(i,j)/radius;
        %
        gr_V(i,j,1)=x_fV(i,j)/radius;
        gr_V(i,j,2)=y_fV(i,j)/radius;
        gr_V(i,j,3)=z_fV(i,j)/radius;
        %
        gr_VI(i,j,1)=x_fVI(i,j)/radius;
        gr_VI(i,j,2)=y_fVI(i,j)/radius;
        gr_VI(i,j,3)=z_fVI(i,j)/radius;
    end
end

%% -- spline data

pmat=zeros(n); % 
pmat=sparse(pmat);
kmat=zeros(n);
kmat=sparse(kmat);
% VALUE OF ETA_CROSS= VALEUR DE ETA_J SUR PANEL II
eta_c=zeros(nn);
for i=1:nn,
    xwk1=sqrt(1+xx(i)^2);
    eta_c(i,1:nn)=atan(xwk1*tan(betacr(i,1:nn)));
end
% ADDRESS OF THE INTERVAL WHERE ETA_C IS LOCATED
jeta_c=zeros(nn);
xwk=zeros(nn-1,1);
eps1=1.e-12;
for j=1:nn,
   for i=1:nn,
       jwk=abs(eta-eta_c(i,j));
       ljwk=(jwk<eps1);
       if sum(ljwk>0),%% eta_c(i,j) coincides with one value of eta
           jeta_c(i,j)=find(ljwk,1);
           jeta_c(i,j)=min(jeta_c(i,j),nn-1);
       else 
         for j1=1:nn-1,
           xwk(j1)=(eta_c(i,j)-eta(j1))*(eta_c(i,j)-eta(j1+1));
         end  
         jeta_c(i,j)=find(xwk<0, 1 );
       end
    end
end
% VALUE OF XI_CROSS= VALEUR DE XI_I SUR PANEL V
xi_c=zeros(nn);
for i=1:nn,
    for j=1:nn,
        xi_c(i,j)=eta_c(j,i);
    end
end
% ADDRESS OF THE INTERVAL WHERE XI_C IS LOCATED
ixi_c=zeros(nn);
xwk=zeros(nn-1,1);
eps1=1.e-12;
for i=1:nn,
    for j=1:nn,
        ixi_c(i,j)=jeta_c(j,i);
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % COEFFICIENTS POUR CONDITION NK
%% % COEFFICIENTS POUR CONDITION AA2
alfasp=1;betasp=2;gamasp=0;
betas=2;gamas=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:n-1,
    pmat(j,j)=4;
    pmat(j,j+1)=1;
    pmat(j,j-1)=1;
end
pmat(1,1)=(4-betasp/alfasp);
pmat(1,2)=(1-gamasp/alfasp);
pmat(n,n-1)=(1-gamasp/alfasp);
pmat(n,n)=(4-betasp/alfasp);
%
%% MATRICE KMAT POUR RELATION HERMITIENNE
for j=2:n-1,
    kmat(j,j+1)=1;
    kmat(j,j-1)=-1;
end
%
kmat(1,1)=-betas/(3*alfasp);
kmat(1,2)=1-gamas/(6*alfasp);
kmat(n,n)=betas/(3*alfasp);
kmat(n,n-1)=-(1-gamas/(6*alfasp));
%%%%%%%%%%%%%%%%
m1=-6+(2*betas+gamas)/alfasp; %% COEFFT LIGNE 1
m2=6-(2*betas+gamas)/alfasp; %% COEFFT LIGNE N-1
%%%%%%%%%%%%%%%%%%%%%%
%% LU FACTORIZATION PF PMAT
% [lmat1,umat1,pmat1,qmat1,rmat1]=lu(pmat);
[lmat1,umat1]=lu(pmat);

%% --- MATRIX FOR DERIVATIVES
% -------------------------------------------------------------------------
pg=zeros(na); 
kg=zeros(na);
for i=2:na-1
    pg(i,i)=4;
    pg(i,i+1)=1;
    pg(i,i-1)=1;
    kg(i,i+1)=1;
    kg(i,i-1)=-1;
end
pg(1,1)=4;pg(1,2)=1;pg(1,na)=1;
pg( na,1)=1;pg(na,na-1)=1;pg(na,na)=4;
kg(1,2)=1;kg(1,na)=-1;
kg(na,1)=1;kg(na,na-1)=-1;
pg=pg./6;
kg=kg./(2*dxi);
% -------------------------------------------------------------------------
if strcmp(scheme,'compact4')==1
    
    for i=2:na-1
        p(i,i)=4;
        p(i,i+1)=1;
        p(i,i-1)=1;
        k(i,i+1)=1;
        k(i,i-1)=-1;
    end
    p(1,1)=4;p(1,2)=1;p(1,na)=1;
    p( na,1)=1;p(na,na-1)=1;p(na,na)=4;
    k(1,2)=1;k(1,na)=-1;
    k(na,1)=1;k(na,na-1)=-1;
    pg=p./6;
    kg=k./(2*deta);

elseif strcmp(scheme,'compact8')==1
    J=diag(ones(na-1,1),-1);
    J(1,end)=1;

    a=25/16;
    b=1/5;
    c=-1/80;
    alpha=3/8;
    k_div=(-a/(2*dxi))*J+(-b/(4*dxi))*J*J+(-c/(6*dxi))*J^3+...
    (a/(2*dxi))*J^(na-1)+(b/(4*dxi))*J^(na-2)+(c/(6*dxi))*J^(na-3);
    kg=sparse(k_div);

    p_div=diag(alpha*ones(na-1,1),1)+diag(alpha*ones(na-1,1),-1)+diag(ones(na,1));
    p_div(1,end)=alpha;
    p_div(end,1)=alpha;
    pg=sparse(p_div);
    
elseif strcmp(scheme,'explicite2')==1
    J=diag(ones(na-1,1),-1);
    J(1,end)=1;

    a=1;
    k=(-a/(2*dxi))*J+(a/(2*dxi))*J^(na-1);
    kg=sparse(k);

    pg=speye(na,na);
    
elseif strcmp(scheme,'explicite4')==1
    J=diag(ones(na-1,1),-1);
    J(1,end)=1;

    a=4/3;
    b=-1/3;
    k=(-a/(2*dxi))*J+(-b/(4*dxi))*J*J+...
    (a/(2*dxi))*J^(na-1)+(b/(4*dxi))*J^(na-2);
    kg=sparse(k);
    pg=speye(na,na);
    
else
    error(['Error in scheme. The spatial scheme : ', scheme ,' is uncorrect'])
end

%% --- Options sur les filtres
% classic filter
if strcmp(opt_ftr,'inf')==1
    f0=1; 
    f1=0;
    f2=0;
    f3=0;
    f4=0;
    f5=0;
    
elseif strcmp(opt_ftr,'redonnet2')==1
    f0=1/2; 
    f1=1/4;
    f2=0;
    f3=0;
    f4=0;
    f5=0;

elseif strcmp(opt_ftr,'redonnet4')==1
    f0=10/16; 
    f1=4/16;
    f2=-1/16;
    f3=0;
    f4=0;
    f5=0;
    
elseif strcmp(opt_ftr,'redonnet6')==1
    f0=44/64; 
    f1=15/64;
    f2=-6/64;
    f3=1/64;
    f4=0;
    f5=0;
    
elseif strcmp(opt_ftr,'redonnet8')==1
    f0=186/256; 
    f1=56/256;
    f2=-28/256;
    f3=8/256;
    f4=-1/256;
    f5=0;

elseif strcmp(opt_ftr,'redonnet10')==1
    f0=772/1024; 
    f1=210/1024;
    f2=-120/1024;
    f3=45/1024;
    f4=-10/1024;
    f5=1/1024;
    
elseif strcmp(opt_ftr,'bogey6')==1
    d0=0.234810479761700;
    d1=-.199250131285813;
    d2=0.120198310245186;
    d3=-.049303775636020;
    d4=0.012396449873964;
    d5=-.001446093078167;
    
    f0=1-d0;
    f1=-d1;
    f2=-d2;
    f3=-d3;
    f4=-d4;
    f5=-d5;
end
%
lig1=[0,1, zeros(1,na-2)];
col1=[zeros(na-1,1);1];
sh1=toeplitz(col1,lig1);
sh1i=inv(sh1);
sh12=sh1^2;
sh1i2=sh1i^2;
sh13=sh12*sh1;
sh1i3=sh1i2*sh1i;
sh14=sh13*sh1;
sh1i4=sh1i3*sh1i;
sh15=sh14*sh1;
sh1i5=sh1i4*sh1i;
ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
ftr=sparse(ftr);

%% INTEGRALE CORRIGEE
% nhs_max=min(127,.125*(nn-1)^2-1); % nombre d'harmoniques à corriger
% [A,err_i] = compute_A_sym(nhs_max); sym=1;
% err=1;
% k=0;
% res = solve_weights( A,err_i,k,err,sym );
% weights=dxi*deta*(dga+res);

weights=dga;