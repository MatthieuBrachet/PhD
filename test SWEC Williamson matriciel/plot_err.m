function []=plot_err(n1,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
global radius;

nn1=n1+2;
mm1=((nn1-1)/2)+1;
%% physical data
radius=6.37122d+06;


%% --- Points on the CS
xi1=linspace(-pi/4, pi/4, nn1); 
dxi1=(pi/2)/(nn1-1);
eta1=linspace(-pi/4, pi/4, nn1); 
deta1=(pi/2)/(nn1-1);
for i=1:nn1,
    xwk=tan(xi1(i));
  for j=1:nn1,
    xx1(i,j)=xwk;
  end
end
for j=1:nn1,
    xwk=tan(eta1(j));
  for i=1:nn1,
    yy1(i,j)=xwk;
  end
end
delta1=zeros(nn1,nn1);
for i=1:nn1,
  for j=1:nn1,
    delta1(i,j)=1+xx1(i,j)^2+yy1(i,j)^2;
  end
end
deltab1=sqrt(delta1); 
dga1=(radius^2)*((1+xx1.^2).*(1+yy1.^2))./(delta1.*deltab1); 
alfa1=zeros(nn1,nn1);
beta1=zeros(nn1,nn1);
for j=1:nn1,
    jbar=j-mm1;
    etabar=jbar*deta1;
    xwk=sqrt(1+tan(etabar)^2);
    for i=1:nn1,
        alfa1(i,j)=atan(tan(xi1(i))/xwk);
    end
end
for i=1:nn1,
    ibar=i-mm1;
    xibar=ibar*dxi1;
    xwk=sqrt(1+tan(xibar)^2);
    for j=1:nn1,
        beta1(i,j)=atan(tan(eta1(j))/xwk);
    end
end
% -------------------------------------------------------------------------
alfacr1=zeros(nn1,nn1);
betacr1=zeros(nn1,nn1);
for j=1:nn1,
    for i=1:nn1,
        betacr1(i,j)=atan(sin(-xi1(i))*tan(eta1(j)));
    end
end
for i=1:nn1
    for j=1:nn1,
        alfacr1(i,j)=acos((cos(betacr1(i,j))/cos(eta1(j)))*sin(-xi1(i)));
    end
end

alfa11=zeros(nn1,nn1); 
for i=1:nn1,
    for j=1:nn1,
        alfa11(i,j)=betacr1(j,i);
    end
end

% -------------------------------------------------------------------------
alfag1=zeros(4*(nn1-1),nn1); 
for j=1:nn1,
  alfag1(1:nn1-1,j)=alfa1(1:nn1-1,j);
  alfag1(nn1:2*(nn1-1),j)=alfacr1(1:nn1-1,j);
  alfag1(2*(nn1-1)+1:3*(nn1-1),j)=alfa1(1:nn1-1,j)+pi;
  alfag1(3*(nn1-1)+1:4*(nn1-1),j)=alfacr1(1:nn1-1,j)+pi;
end
betag1=alfag1'; 
% - FACE I -
for i=1:nn1,  
    for j=1:nn1,
        jbar=j-mm1;
        etabar=jbar*deta1;
        al1=alfa1(i,j);
        x_I(i,j)=cos(al1)*cos(etabar);
        y_I(i,j)=sin(al1);
        z_I(i,j)=cos(al1)*sin(etabar);
    end
end
x_I=radius*x_I;
y_I=radius*y_I;
z_I=radius*z_I;
% - FACE II -
for i=1:nn1, 
    for j=1:nn1,
        jbar=j-mm1;
        etabar=jbar*deta1;
        al2=alfa1(i,j);
        x_II(i,j)=-sin(al2);
        y_II(i,j)=cos(al2)*cos(etabar);
        z_II(i,j)=cos(al2)*sin(etabar);
    end
end
x_II=radius*x_II;
y_II=radius*y_II;
z_II=radius*z_II;
% - FACE III -
for i=1:nn1,  
    for j=1:nn1,
        jbar=j-mm1;
        etabar=jbar*deta1;
        al3=alfa1(i,j);
        x_III(i,j)=-cos(al3)*cos(etabar);
        y_III(i,j)=-sin(al3);
        z_III(i,j)=cos(al3)*sin(etabar);
    end
end
x_III=radius*x_III;
y_III=radius*y_III;
z_III=radius*z_III;
% - FACE IV -
for i=1:nn1,  
    for j=1:nn1,
        jbar=j-mm1;
        etabar=jbar*deta1;
        al4=alfa1(i,j);
        x_IV(i,j)=sin(al4);
        y_IV(i,j)=-cos(al4)*cos(etabar);
        z_IV(i,j)=cos(al4)*sin(etabar);
    end
end
x_IV=radius*x_IV;
y_IV=radius*y_IV;
z_IV=radius*z_IV;



%% *** PLOT ***************************************************************

%% --- axis
lmin=-pi;lmax=pi;temin=-pi/2;temax=pi/2;
umin=min(min([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));
umax=max(max([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));

%% -- FACE I
[lam_I,the_I,xwk]=cart2sph(x_I,y_I,z_I);
surf(lam_I,the_I,funfIe);hold on;axis([lmin lmax temin temax umin umax]);

%% -- FACE II
[lam_II,the_II,xwk]=cart2sph(x_II,y_II,z_II);
surf(lam_II,the_II,funfIIe);hold on;axis([lmin lmax temin temax umin umax]);

%% -- FACE III
[lam_III,the_III,xwk]=cart2sph(x_III,y_III,z_III);
lam_IIIa=lam_III((nn1+1)/2:nn1,1:nn1);
the_IIIa=the_III((nn1+1)/2:nn1,1:nn1);
funfIIIea=funfIIIe((nn1+1)/2:nn1,1:nn1);
surf(lam_IIIa,the_IIIa,funfIIIea);hold on;axis([lmin lmax temin temax umin umax]);
lam_IIIb=lam_III(1:(nn1+1)/2,1:nn1);
lam_IIIb=lam_IIIb+2*pi*(lam_IIIb<-eps);
the_IIIb=the_III(1:(nn1+1)/2,1:nn1);
funfIIIeb=funfIIIe(1:(nn1+1)/2,1:nn1);
surf(lam_IIIb,the_IIIb,funfIIIeb);hold on;axis([lmin lmax temin temax umin umax]);

%% -- FACE IV
[lam_IV,the_IV,xwk]=cart2sph(x_IV,y_IV,z_IV);
surf(lam_IV,the_IV,funfIVe);hold on;axis([lmin lmax temin temax umin umax]);

%% -- FACE V
% QUADRANT I
funfVsI=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2 
    alfaspline(1:nn1)=alfa1(1:nn1,j);
    funspline(1:nn1)=funfVe(1:nn1,j);
    ppspline=spline(alfaspline,funspline);
    funfVsI(1:nn1,j)=ppval(ppspline,alfa11(1:nn1,j));
end
lam_V_qI=zeros(nn1,(nn1+1)/2);
the_V_qI=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    etawk=eta1(j);
    for i=1:nn1 
         alfawk=alfa11(i,j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qI(i,j),the_V_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qI(1:nn1,(nn1+1)/2)=lam_V_qI(1:nn1,1);
surf(lam_V_qI,the_V_qI,funfVsI); hold on; axis([lmin lmax temin temax umin umax]);

% QUADRANT II
funfVsII=zeros((nn1+1)/2,nn1);
for i=(nn1+1)/2:-1:1
    betaspline(1:nn1)=beta1(i+(nn1-1)/2,1:nn1);
    funspline(1:nn1)=funfVe(i+(nn1-1)/2,1:nn1);
    ppspline=spline(betaspline,funspline);
    funfVsII(i,1:nn1)=ppval(ppspline,betacr1(i+(nn1-1)/2,1:nn1));
end
lam_V_qII=zeros(nn1,(nn1+1)/2);
the_V_qII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2
    xiwk=xi1(nn1-j+1);
    for i=1:nn1
         betawk = -betacr1(j,i);
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qII(i,j),the_V_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qII(1:nn1,(nn1+1)/2)=lam_V_qII(1:nn1,1);
funfVsII1=funfVsII((nn1+1)/2:-1:1,:)';
surf(lam_V_qII,the_V_qII,funfVsII1);hold on;axis([lmin lmax temin temax umin umax]);

% QUADRANT III
funfVsIII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    alfaspline(1:nn1)=alfa1(1:nn1,j);
    funspline(1:nn1)=funfVe(1:nn1,nn1+1-j);
    ppspline=spline(alfaspline,funspline);
    funfVsIII(1:nn1,j)=ppval(ppspline,alfa11(nn1:-1:1,nn1+1-j));
end
lam_V_qIII=zeros(nn1,(nn1+1)/2);
the_V_qIII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    etawk=eta1(j);
    for i=1:nn1
         alfawk=alfa11(i,nn1+1-j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qIII(i,j),the_V_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIII(1:nn1,(nn1+1)/2)=lam_V_qIII(1:nn1,1);
lam_V_a=lam_V_qIII((nn1+1)/2:nn1,:)+pi;
the_V_a=the_V_qIII((nn1+1)/2:nn1,:);
funfVa=funfVsIII((nn1+1)/2:nn1,:);
surf(lam_V_a,the_V_a,funfVa);hold on;axis([lmin lmax temin temax umin umax]);
lam_V_b=lam_V_qIII(1:(nn1+1)/2,:)-pi;
the_V_b=the_V_qIII(1:(nn1+1)/2,:);
funfVb=funfVsIII(1:(nn1+1)/2,:);
surf(lam_V_b,the_V_b,funfVb);hold on;axis([lmin lmax temin temax umin umax]);

% QUADRANT IV
funfVsIV=zeros((nn1+1)/2,nn1); % 
for i=1:(nn1+1)/2,
    betaspline(1:nn1)=beta1(i,1:nn1);
    funspline(1:nn1)=funfVe(i,1:nn1);
    ppspline=spline(betaspline,funspline);
     funfVsIV(i,1:nn1)=ppval(ppspline,betacr1(i,nn1:-1:1));
end
lam_V_qIV=zeros(nn1,(nn1+1)/2);
the_V_qIV=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    xiwk=xi1(j);
    for i=1:nn1 
         betawk = betacr1(j,nn1-i+1);
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qIV(i,j),the_V_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIV(1:nn1,(nn1+1)/2)=lam_V_qIV(1:nn1,1);
funfVsIV1=funfVsIV(1:(nn1+1)/2,:)';
lam_V_qIV1=lam_V_qIV+(2*pi)*(lam_V_qIV<0);
surf(lam_V_qIV1-2*pi,the_V_qIV,funfVsIV1);hold on;axis([lmin lmax temin temax umin umax]);

%% -- PANEL VI
% QUADRANT I
funfVIsI=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2
    alfaspline(1:nn1)=alfa1(1:nn1,j+(nn1+1)/2-1);
    funspline(1:nn1)=funfVIe(1:nn1,j+(nn1+1)/2-1);
    ppspline=spline(alfaspline,funspline);
    funfVIsI(1:nn1,j)=ppval(ppspline,alfa11(1:nn1,j+(nn1+1)/2-1));
end
lam_VI_qI=zeros(nn1,(nn1+1)/2);
the_VI_qI=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    etawk=eta1(j+(nn1+1)/2-1); 
    for i=1:nn1
         alfawk=alfa11(i,(nn1+1)/2-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qI(i,j),the_VI_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qI(1:nn1,1)=lam_VI_qI(1:nn1,(nn1+1)/2);
surf(lam_VI_qI,the_VI_qI,funfVIsI(nn1:-1:1,:));hold on;axis([lmin lmax temin temax umin umax]);

% QUADRANT II
funfVIsII=zeros((nn1+1)/2,nn1);
for i=(nn1+1)/2:-1:1,
    betaspline(1:nn1)=beta1(i+(nn1-1)/2,1:nn1);
    funspline(1:nn1)=funfVIe(i+(nn1-1)/2,1:nn1);
    ppspline=spline(betaspline,funspline);
    funfVIsII(i,1:nn1)=ppval(ppspline,betacr1(i+(nn1-1)/2,1:nn1));
end
lam_VI_qII=zeros(nn1,(nn1+1)/2);
the_VI_qII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2,
    xiwk=xi1(nn1-j+1);
    for i=1:nn1 
         betawk = betacr1(j,i); 
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qII(i,j),the_VI_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qII(1:nn1,(nn1+1)/2)=lam_VI_qII(1:nn1,1);
funfVIsII1=funfVIsII((nn1+1)/2:-1:1,nn1:-1:1)';
surf(lam_VI_qII,the_VI_qII,funfVIsII1);hold on;axis([lmin lmax temin temax umin umax]);

% QUADRANT III
funfVIsIII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2
    alfaspline(1:nn1)=alfa1(1:nn1,j);
    funspline(1:nn1)=funfVIe(1:nn1,j);
    ppspline=spline(alfaspline,funspline);
    funfVIsIII(1:nn1,j)=ppval(ppspline,alfa11(nn1:-1:1,j));
end
lam_VI_qIII=zeros(nn1,(nn1+1)/2);
the_VI_qIII=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2
    etawk=eta1(nn1-j+1);
    for i=1:nn1
         alfawk=alfa11(nn1-i+1,nn1-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qIII(i,j),the_VI_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIII(1:nn1,(nn1+1)/2)=lam_VI_qIII(1:nn1,1);
lam_VI_a=lam_VI_qIII(1:(nn1+1)/2,:)+pi;
the_VI_a=the_VI_qIII(1:(nn1+1)/2,:);
funfVIa=funfVIsIII(1:(nn1+1)/2,:);
surf(lam_VI_a,the_VI_a,funfVIa);hold on;axis([lmin lmax temin temax umin umax]);
lam_VI_b=lam_VI_qI((nn1+1)/2:nn1,:)-pi;
the_VI_b=the_VI_qIII((nn1+1)/2:nn1,:);
funfVIb=funfVIsIII((nn1+1)/2:nn1,:);
surf(lam_VI_b,the_VI_b,funfVIb);hold on;axis([lmin lmax temin temax umin umax]);

% QUADRANT IV
funfVIsIV=zeros((nn1+1)/2,nn1);
for i=1:(nn1+1)/2,
    betaspline(1:nn1)=beta1(i,1:nn1);
    funspline(1:nn1)=funfVIe(i,1:nn1);
    ppspline=spline(betaspline,funspline);
    funfVIsIV(i,1:nn1)=ppval(ppspline,betacr1(i,nn1:-1:1));
end
lam_VI_qIV=zeros(nn1,(nn1+1)/2);
the_VI_qIV=zeros(nn1,(nn1+1)/2);
for j=1:(nn1+1)/2
    xiwk=xi1((nn1+1)/2-j+1);
    for i=1:nn1
         betawk = betacr1((nn1+1)/2-j+1,nn1-i+1);
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qIV(i,j),the_VI_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIV(1:nn1,1)=lam_VI_qIV(1:nn1,(nn1+1)/2);
funfVIsIV=funfVIsIV((nn1+1)/2:-1:1,1:nn1)';
surf(lam_VI_qIV,the_VI_qIV,funfVIsIV);hold on;axis([lmin lmax temin temax umin umax]);


%% -- TRACE DES FRONTIERES DES PATCHS DE LA CUBED SPHERE GRID
plot3(lam_I(1,:),the_I(1,:),funfIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(nn1,:),the_I(nn1,:),funfIe(nn1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(:,1),the_I(:,1),funfIe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(:,nn1),the_I(:,nn1),funfIe(:,nn1)+eps,'k','LineWidth',1.25); hold on;

plot3(lam_II(1,:),the_II(1,:),funfIIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(nn1,:),the_II(nn1,:),funfIIe(nn1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,1),the_II(:,1),funfIIe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,nn1),the_II(:,nn1),funfIIe(:,nn1)+eps,'k','LineWidth',1.25); hold on;

lam_III=lam_III+2*pi*(lam_III<-eps);
plot3(lam_III(1,:),the_III(1,:),funfIIIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(nn1,:),the_III(nn1,:),funfIIIe(nn1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(1:(nn1+1)/2,1),the_III(1:(nn1+1)/2,1),funfIIIe(1:(nn1+1)/2,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III((nn1+1)/2:nn1,1)-2*pi,the_III((nn1+1)/2:nn1,1),funfIIIe((nn1+1)/2:nn1,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(1:(nn1+1)/2,nn1),the_III(1:(nn1+1)/2,nn1),funfIIIe(1:(nn1+1)/2,nn1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III((nn1+1)/2:nn1,nn1)-2*pi,the_III((nn1+1)/2:nn1,nn1),funfIIIe((nn1+1)/2:nn1,nn1)+eps,'k','LineWidth',1.25); hold on;

plot3(lam_IV(1,:),the_IV(1,:),funfIVe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(nn1,:),the_IV(nn1,:),funfIVe(nn1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,1),the_IV(:,1),funfIVe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,nn1),the_IV(:,nn1),funfIVe(:,nn1)+eps,'k','LineWidth',1.25); hold on;

umax=max(max([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));
load coast
sz=size(long);
plot3(long*pi/180,lat*pi/180,umax*ones(sz),'k-')

shading interp;
set(gca, 'CLim', [umin, umax]);

view(2);






end