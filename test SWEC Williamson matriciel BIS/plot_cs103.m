function []=plot_cs103(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,v)
% classic plot en projection stéréographique en centrant le 0 sur
% la figure
% lignes de niveaux suivant v
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global xi eta
global alfa beta;
global betacr;
global alfa1;

%% --- axis
lmin=-pi;lmax=pi;temin=-pi/2;temax=pi/2;
umin=min(min([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));
umax=max(max([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));

%% -- FACE I
[lam_I,the_I,xwk]=cart2sph(x_fI,y_fI,z_fI);
contourf(lam_I,the_I,funfIe.*(funfIe>0),v,'k-');hold on;axis([lmin lmax temin temax]);
contourf(lam_I,the_I,funfIe.*(funfIe<0),v,'k--');hold on;

%% -- FACE II
[lam_II,the_II,xwk]=cart2sph(x_fII,y_fII,z_fII);
contourf(lam_II,the_II,funfIIe.*(funfIIe>0),v,'k-');hold on;
contourf(lam_II,the_II,funfIIe.*(funfIIe<0),v,'k--');hold on;

%% -- FACE III
[lam_III,the_III,xwk]=cart2sph(x_fIII,y_fIII,z_fIII);
lam_IIIa=lam_III((nn+1)/2:nn,1:nn);
the_IIIa=the_III((nn+1)/2:nn,1:nn);
funfIIIea=funfIIIe((nn+1)/2:nn,1:nn);
contourf(lam_IIIa,the_IIIa,funfIIIea.*(funfIIIea>0),v,'k-');hold on;
contourf(lam_IIIa,the_IIIa,funfIIIea.*(funfIIIea<0),v,'k--');hold on;

lam_IIIb=lam_III(1:(nn+1)/2,1:nn);
lam_IIIb=lam_IIIb+2*pi*(lam_IIIb<-eps);
the_IIIb=the_III(1:(nn+1)/2,1:nn);
funfIIIeb=funfIIIe(1:(nn+1)/2,1:nn);
contourf(lam_IIIb,the_IIIb,funfIIIeb.*(funfIIIeb>0),v,'k-');hold on;
contourf(lam_IIIb,the_IIIb,funfIIIeb.*(funfIIIeb<0),v,'k--');hold on;

%% -- FACE IV
[lam_IV,the_IV,xwk]=cart2sph(x_fIV,y_fIV,z_fIV);
contourf(lam_IV,the_IV,funfIVe.*(funfIVe>0),v,'k-');hold on;
contourf(lam_IV,the_IV,funfIVe.*(funfIVe<0),v,'k--');hold on;

%% -- FACE V
% QUADRANT I
funfVsI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVe(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funfVsI(1:nn,j)=ppval(ppspline,alfa1(1:nn,j));
end
lam_V_qI=zeros(nn,(nn+1)/2);
the_V_qI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    etawk=eta(j);
    for i=1:nn 
         alfawk=alfa1(i,j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qI(i,j),the_V_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qI(1:nn,(nn+1)/2)=lam_V_qI(1:nn,1);
contourf(lam_V_qI,the_V_qI,funfVsI.*(funfVsI>0),v,'k-'); hold on;
contourf(lam_V_qI,the_V_qI,funfVsI.*(funfVsI<0),v,'k--'); hold on;

% QUADRANT II
funfVsII=zeros((nn+1)/2,nn);
for i=(nn+1)/2:-1:1
    betaspline(1:nn)=beta(i+(nn-1)/2,1:nn);
    funspline(1:nn)=funfVe(i+(nn-1)/2,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVsII(i,1:nn)=ppval(ppspline,betacr(i+(nn-1)/2,1:nn));
end
lam_V_qII=zeros(nn,(nn+1)/2);
the_V_qII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
    xiwk=xi(nn-j+1);
    for i=1:nn
         betawk = -betacr(j,i);
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qII(i,j),the_V_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qII(1:nn,(nn+1)/2)=lam_V_qII(1:nn,1);
funfVsII1=funfVsII((nn+1)/2:-1:1,:)';
contourf(lam_V_qII,the_V_qII,funfVsII1.*(funfVsII1>0),v,'k-');hold on;
contourf(lam_V_qII,the_V_qII,funfVsII1.*(funfVsII1<0),v,'k--');hold on;

% QUADRANT III
funfVsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVe(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funfVsIII(1:nn,j)=ppval(ppspline,alfa1(nn:-1:1,nn+1-j));
end
lam_V_qIII=zeros(nn,(nn+1)/2);
the_V_qIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    etawk=eta(j);
    for i=1:nn
         alfawk=alfa1(i,nn+1-j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qIII(i,j),the_V_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIII(1:nn,(nn+1)/2)=lam_V_qIII(1:nn,1);
lam_V_a=lam_V_qIII((nn+1)/2:nn,:)+pi;
the_V_a=the_V_qIII((nn+1)/2:nn,:);
funfVa=funfVsIII((nn+1)/2:nn,:);
contourf(lam_V_a,the_V_a,funfVa.*(funfVa>0),v,'k-');hold on;
contourf(lam_V_a,the_V_a,funfVa.*(funfVa<0),v,'k--');hold on;

lam_V_b=lam_V_qIII(1:(nn+1)/2,:)-pi;
the_V_b=the_V_qIII(1:(nn+1)/2,:);
funfVb=funfVsIII(1:(nn+1)/2,:);
contourf(lam_V_b,the_V_b,funfVb.*(funfVb>0),v,'k-');hold on;
contourf(lam_V_b,the_V_b,funfVb.*(funfVb<0),v,'k--');hold on;

% QUADRANT IV
funfVsIV=zeros((nn+1)/2,nn); 
for i=1:(nn+1)/2,
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVe(i,1:nn);
    ppspline=spline(betaspline,funspline);
     funfVsIV(i,1:nn)=ppval(ppspline,betacr(i,nn:-1:1));
end
lam_V_qIV=zeros(nn,(nn+1)/2);
the_V_qIV=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    xiwk=xi(j);
    for i=1:nn 
         betawk = betacr(j,nn-i+1);
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qIV(i,j),the_V_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIV(1:nn,(nn+1)/2)=lam_V_qIV(1:nn,1);
funfVsIV1=funfVsIV(1:(nn+1)/2,:)';
lam_V_qIV1=lam_V_qIV+(2*pi)*(lam_V_qIV<0);
contourf(lam_V_qIV1-2*pi,the_V_qIV,funfVsIV1.*(funfVsIV1>0),v,'k-');hold on;
contourf(lam_V_qIV1-2*pi,the_V_qIV,funfVsIV1.*(funfVsIV1<0),v,'k--');hold on;

%% -- PANEL VI
% QUADRANT I
funfVIsI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
    alfaspline(1:nn)=alfa(1:nn,j+(nn+1)/2-1);
    funspline(1:nn)=funfVIe(1:nn,j+(nn+1)/2-1);
    ppspline=spline(alfaspline,funspline);
    funfVIsI(1:nn,j)=ppval(ppspline,alfa1(1:nn,j+(nn+1)/2-1));
end
lam_VI_qI=zeros(nn,(nn+1)/2);
the_VI_qI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    etawk=eta(j+(nn+1)/2-1); 
    for i=1:nn
         alfawk=alfa1(i,(nn+1)/2-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qI(i,j),the_VI_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qI(1:nn,1)=lam_VI_qI(1:nn,(nn+1)/2);
funf=funfVIsI(nn:-1:1,:);
contourf(lam_VI_qI,the_VI_qI,funf.*(funf>0),v,'k-');hold on;
contourf(lam_VI_qI,the_VI_qI,funf.*(funf<0),v,'k--');hold on;

% QUADRANT II
funfVIsII=zeros((nn+1)/2,nn);
for i=(nn+1)/2:-1:1,
    betaspline(1:nn)=beta(i+(nn-1)/2,1:nn);
    funspline(1:nn)=funfVIe(i+(nn-1)/2,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVIsII(i,1:nn)=ppval(ppspline,betacr(i+(nn-1)/2,1:nn));
end
lam_VI_qII=zeros(nn,(nn+1)/2);
the_VI_qII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    xiwk=xi(nn-j+1);
    for i=1:nn 
         betawk = betacr(j,i); 
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qII(i,j),the_VI_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qII(1:nn,(nn+1)/2)=lam_VI_qII(1:nn,1);
funfVIsII1=funfVIsII((nn+1)/2:-1:1,nn:-1:1)';
contourf(lam_VI_qII,the_VI_qII,funfVIsII1.*(funfVIsII1>0),v,'k-');hold on;
contourf(lam_VI_qII,the_VI_qII,funfVIsII1.*(funfVIsII1<0),v,'k--');hold on;

% QUADRANT III
funfVIsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVIe(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funfVIsIII(1:nn,j)=ppval(ppspline,alfa1(nn:-1:1,j));
end
lam_VI_qIII=zeros(nn,(nn+1)/2);
the_VI_qIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
    etawk=eta(nn-j+1);
    for i=1:nn
         alfawk=alfa1(nn-i+1,nn-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qIII(i,j),the_VI_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIII(1:nn,(nn+1)/2)=lam_VI_qIII(1:nn,1);
lam_VI_a=lam_VI_qIII(1:(nn+1)/2,:)+pi;
the_VI_a=the_VI_qIII(1:(nn+1)/2,:);
funfVIa=funfVIsIII(1:(nn+1)/2,:);
contourf(lam_VI_a,the_VI_a,funfVIa.*(funfVIa>0),v,'k-');hold on;
contourf(lam_VI_a,the_VI_a,funfVIa.*(funfVIa<0),v,'k--');hold on;

lam_VI_b=lam_VI_qI((nn+1)/2:nn,:)-pi;
the_VI_b=the_VI_qIII((nn+1)/2:nn,:);
funfVIb=funfVIsIII((nn+1)/2:nn,:);
contourf(lam_VI_b,the_VI_b,funfVIb.*(funfVIb>0),v,'k-');hold on;
contourf(lam_VI_b,the_VI_b,funfVIb.*(funfVIb<0),v,'k--');hold on;

% QUADRANT IV
funfVIsIV=zeros((nn+1)/2,nn);
for i=1:(nn+1)/2,
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVIe(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVIsIV(i,1:nn)=ppval(ppspline,betacr(i,nn:-1:1));
end
lam_VI_qIV=zeros(nn,(nn+1)/2);
the_VI_qIV=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
    xiwk=xi((nn+1)/2-j+1);
    for i=1:nn 
         betawk = betacr((nn+1)/2-j+1,nn-i+1);
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qIV(i,j),the_VI_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIV(1:nn,1)=lam_VI_qIV(1:nn,(nn+1)/2);
funfVIsIV=funfVIsIV((nn+1)/2:-1:1,1:nn)';
contourf(lam_VI_qIV,the_VI_qIV,funfVIsIV.*(funfVIsIV>0),v,'k-');hold on;
contourf(lam_VI_qIV,the_VI_qIV,funfVIsIV.*(funfVIsIV<0),v,'k--');hold on;


%% -- TRACE DES FRONTIERES DES PATCHS DE LA CUBED SPHERE GRID
lam_Ia=lam_I((nn+1)/2:nn,1:nn);
the_Ia=the_I((nn+1)/2:nn,1:nn);

lam_Ib=lam_I(1:(nn+1)/2,1:nn);
the_Ib=the_I(1:(nn+1)/2,1:nn);

funfIea=funfIe((nn+1)/2:nn,1:nn);
funfIeb=funfIe(1:(nn+1)/2,1:nn);

plot3(lam_I(1,:),the_I(1,:),funfIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(nn,:),the_I(nn,:),funfIe(nn,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(:,1),the_I(:,1),funfIe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_I(:,nn),the_I(:,nn),funfIe(:,nn)+eps,'k','LineWidth',1.25); hold on;

plot3(lam_II(1,:),the_II(1,:),funfIIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(nn,:),the_II(nn,:),funfIIe(nn,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,1),the_II(:,1),funfIIe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,nn),the_II(:,nn),funfIIe(:,nn)+eps,'k','LineWidth',1.25); hold on;

lam_III=lam_III+2*pi*(lam_III<-eps);
plot3(lam_III(1,:),the_III(1,:),funfIIIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(nn,:),the_III(nn,:),funfIIIe(nn,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(1:(nn+1)/2,1),the_III(1:(nn+1)/2,1),funfIIIe(1:(nn+1)/2,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III((nn+1)/2:nn,1)-2*pi,the_III((nn+1)/2:nn,1),funfIIIe((nn+1)/2:nn,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III(1:(nn+1)/2,nn),the_III(1:(nn+1)/2,nn),funfIIIe(1:(nn+1)/2,nn)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_III((nn+1)/2:nn,nn)-2*pi,the_III((nn+1)/2:nn,nn),funfIIIe((nn+1)/2:nn,nn)+eps,'k','LineWidth',1.25); hold on;

plot3(lam_IV(1,:),the_IV(1,:),funfIVe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(nn,:),the_IV(nn,:),funfIVe(nn,:)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,1),the_IV(:,1),funfIVe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,nn),the_IV(:,nn),funfIVe(:,nn)+eps,'k','LineWidth',1.25); hold on;


shading interp;
set(gca, 'CLim', [umin, umax]);

view(2);