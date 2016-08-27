% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
global n nn;
global mm na nb;
global radius;
global xi eta dxi deta xx yy delta deltab dga;
global alfa beta;
global alfacr betacr;
global alfa1;
global alfag betag;
global p k kxi keta k_div p_div;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global opt_ftr ftr opt_ftr_div ftr_div;
global omega hp gp u0 h0

%% physical data
radius=6.37122d+06;
omega=7.292d-05;
hp=100;
gp=9.80616;
u0=80;
h0=100;


nn=n+2;
% -----------------------------------------
%global mm na nb;
mm=((nn-1)/2)+1;
na=4*(nn-1);
nb=na;
% ----------------------------------------
%global radius;
%radius=6.37122d+06; % rayon terrestre
% -----------------------------------------
%global xi eta dxi deta xx yy delta deltab;
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
deltab=sqrt(delta); % DELTAB=DELTA DE ULLRICH = SQRT(1+X^2+Y^2).
dga=zeros(nn,nn);
dga=(radius^2)*((1+xx.^2).*(1+yy.^2))./(delta.*deltab); % ELEMENT AREA
% ---------------------------------------
%global alfa beta;
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
% -----------------------------------------
%global alfacr betacr;
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
%global alfa1;
alfa1=zeros(nn,nn); % VARIABLE DE TRAVAIL= ANGLE ALPHA DE TYPE BETACROSS
% LIEN ENTRE ALFA1 ET BETACROSS, CE CALCUL EST A REFAIRE POUR ETRE SUR. 
for i=1:nn,
    for j=1:nn,
        alfa1(i,j)=betacr(j,i);
    end
end
% ----------------------------------------
% CALCUL DES ANGLES GLOBAUX DE RESEAUX: 2 TABLEAUX SEULEMENT ALFA_G ET BETA_G
% ALFA_G: ABSCISSES CURVILIGNES LE LONG DE Ia, IIa, Va
% BETA_G: ABSCISSES CURVILIGNES LE LONG DE Ib, IIb, Vb
%global alfag betag;
alfag=zeros(4*(nn-1),nn); 
for j=1:nn,
  alfag(1:nn-1,j)=alfa(1:nn-1,j);
  alfag(nn:2*(nn-1),j)=alfacr(1:nn-1,j);
  alfag(2*(nn-1)+1:3*(nn-1),j)=alfa(1:nn-1,j)+pi;
  alfag(3*(nn-1)+1:4*(nn-1),j)=alfacr(1:nn-1,j)+pi;
end
betag=alfag'; % BETAG=TRANSPOSEE DE ALFAG
%% --------------------------------------------------
%global p k;
p=zeros(na); % 
p=sparse(p);
k=zeros(na);
k=sparse(k);
%
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

p=p/6;
kxi=k./(2*dxi);
keta=k./(2*deta);

J=diag(ones(na-1,1),-1);
J(1,end)=1;
a=25/16;
b=1/5;
c=-1/80;
alpha=3/8;
k_div=(-a/(2*dxi))*J+(-b/(4*dxi))*J*J+(-c/(6*dxi))*J^3+...
    (a/(2*dxi))*J^(na-1)+(b/(4*dxi))*J^(na-2)+(c/(6*dxi))*J^(na-3);
k_div=sparse(k_div);

p_div=diag(alpha*ones(na-1,1),1)+diag(alpha*ones(na-1,1),-1)+diag(ones(na,1));
p_div(1,end)=alpha;
p_div(end,1)=alpha;
p_div=sparse(p_div);

%p_div=p;
%k_div=kxi;



%% ----------------------------------------------------------------
% % CARTESIAN COORDINATES OF THE POINTS OF THE 6 FACES.
% global x_fI y_fI z_fI;
% global x_fII y_fII z_fII;
% global x_fIII y_fIII z_fIII;
% global x_fIV y_fIV z_fIV;
% global x_fV y_fV z_fV;
% global x_fVI y_fVI z_fVI;
% -------------------------------------------------------------------
% MODIFICATION CALCUL DES COORDONNEES CARTESIENNES: PLUS CLAIR.
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
% -----------------------------------------------------------------
% DATA OF THE CONTRAVARIANT VECTORS XI/ETA FOR EACH FACE
% - FACE I -
% global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
% global geta_I geta_II geta_III geta_IV geta_V geta_VI;
gxi_I=zeros(nn,nn,3);geta_I=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
%
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
% - FACE II -
 gxi_II=zeros(nn,nn,3);geta_II=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
% - FACE III -
 gxi_III=zeros(nn,nn,3);geta_III=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
% - FACE IV -
 gxi_IV=zeros(nn,nn,3);geta_IV=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
% - FACE V -
 gxi_V=zeros(nn,nn,3);geta_V=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
% - FACE VI -
 gxi_VI=zeros(nn,nn,3);geta_VI=zeros(nn,nn,3); % CONTRAVARIANT VECTORS XI/ETA
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
% p1=zeros(n); % 
% p1=sparse(p1);
% k1=zeros(n);
% k1=sparse(k1);
% %
% for i=2:n-1
%     p1(i,i)=4;
%     p1(i,i+1)=1;
%     p1(i,i-1)=1;
%     k1(i,i+1)=1;
%     k1(i,i-1)=-1;
% end
% p1(1,1)=4;p1(1,2)=1;
% p1(n,n-1)=1;p1(n,n)=4;
% k1(1,2)=1;
% k1(n,n-1)=-1;


%% Options sur les filtres
[ ftr ] = filtre( na , opt_ftr );
[ ftr_div ] = filtre( na , opt_ftr_div );