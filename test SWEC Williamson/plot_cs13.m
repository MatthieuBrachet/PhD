function []=plot_cs13(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
% projection stereographique et tracé des contours.
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global xi eta
global alfa beta;
global betacr;
global alfa1;

MI=max(max(funfIe));
MII=max(max(funfIIe));
MIII=max(max(funfIIIe));
MIV=max(max(funfIVe));
MV=max(max(funfVe));
MVI=max(max(funfVIe));
MAX=max([MI,MII,MIII,MIV,MV,MVI]);
MI=min(min(funfIe));
MII=min(min(funfIIe));
MIII=min(min(funfIIIe));
MIV=min(min(funfIVe));
MV=min(min(funfVe));
MVI=min(min(funfVIe));
MIN=min([MI,MII,MIII,MIV,MV,MVI]);
v=linspace(MAX,MIN,100);
dat=(MAX+MIN)/2;


 % axis
 % DESSIN DE TYPE CYLINDRIQUE 
 [lam_I,the_I,rwk]=cart2sph(x_fI,y_fI,z_fI);
 [lam_II,the_II,rwk]=cart2sph(x_fII,y_fII,z_fII);
 [lam_III,the_III,rwk]=cart2sph(x_fIII,y_fIII,z_fIII);
 [lam_IV,the_IV,rwk]=cart2sph(x_fIV,y_fIV,z_fIV);
 [lam_V,the_V,rwk]=cart2sph(x_fV,y_fV,z_fV);
 [lam_VI,the_VI,rwk]=cart2sph(x_fVI,y_fVI,z_fVI);
 
 lam_I=lam_I+2*pi*(lam_I<-eps); 
 lam_II=lam_II+2*pi*(lam_II<-eps);
 lam_III=lam_III+2*pi*(lam_III<-eps);
 lam_IV=lam_IV+2*pi*(lam_IV<-eps);
 lam_V=lam_V+2*pi*(lam_V<-eps);
 lam_VI=lam_VI+2*pi*(lam_VI<-eps);
 
 the_I=the_I*180/pi;
 the_II=the_II*180/pi;
 the_III=the_III*180/pi;
 the_IV=the_IV*180/pi;
 the_V=the_V*180/pi;
 the_VI=the_VI*180/pi;
 
 %%
lam_Ia=lam_I((nn+1)/2:nn,1:nn);
the_Ia=the_I((nn+1)/2:nn,1:nn);

lam_Ib=lam_I(1:(nn+1)/2,1:nn);
lam_Ib((nn+1)/2,:)=lam_Ib((nn+1)/2,:)+2*pi;
the_Ib=the_I(1:(nn+1)/2,1:nn);

funfIea=funfIe((nn+1)/2:nn,1:nn);
funfIeb=funfIe(1:(nn+1)/2,1:nn);

%% ************************************************************************
% face I-II-III-IV
%% ************************************************************************
lmin=0;lmax=360;temin=-90;temax=90;
umin=min(min([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));
umax=max(max([funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe]));
axis([lmin lmax temin temax umin umax]);

M=funfIea;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_Ia*180/pi,the_Ia,Mneg,v,'k--');
contour(lam_Ia*180/pi,the_Ia,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

M=funfIIe;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_II*180/pi,the_II,Mneg,v,'k--');
contour(lam_II*180/pi,the_II,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

M=funfIIIe;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_III*180/pi,the_III,Mneg,v,'k--');
contour(lam_III*180/pi,the_III,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]); 

M=funfIVe;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_IV*180/pi,the_IV,Mneg,v,'k--');
contour(lam_IV*180/pi,the_IV,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

M=funfIeb;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_Ib*180/pi,the_Ib,Mneg,v,'k--');
contour(lam_Ib*180/pi,the_Ib,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

%% ************************************************************************
% face V
%% ************************************************************************

%% QUADRANT I OF FACE V
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
lam_V_a=lam_V_qI((nn+1)/2:nn,:);
the_V_a=the_V_qI((nn+1)/2:nn,:);
funfVa=funfVsI((nn+1)/2:nn,:);


M=funfVa;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_V_a*180/pi,the_V_a*180/pi,Mneg,v,'k--');
contour(lam_V_a*180/pi,the_V_a*180/pi,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

lam_V_b=lam_V_qI(1:(nn+1)/2,:)+2*pi;
the_V_b=the_V_qI(1:(nn+1)/2,:);
funfVb=funfVsI(1:(nn+1)/2,:);

M=funfVb;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_V_b*180/pi,the_V_b*180/pi,Mneg,v,'k--');
contour(lam_V_b*180/pi,the_V_b*180/pi,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);


%% QUADRANT III OF FACE V ; REVOIR ADRESSAGE: PAS CLAIR

funfVsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2
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

M=funfVsIII;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour((lam_V_qIII+pi)*180/pi,the_V_qIII*180/pi,Mneg,v,'k--');
contour((lam_V_qIII+pi)*180/pi,the_V_qIII*180/pi,Mpos,v,'k-');
hold on;axis([lmin lmax temin temax umin umax]);

%%     QUADRANT II DE FACE V

funfVsII=zeros((nn+1)/2,nn); 
for i=(nn+1)/2:-1:1,
    betaspline(1:nn)=beta(i+(nn-1)/2,1:nn);
    funspline(1:nn)=funfVe(i+(nn-1)/2,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVsII(i,1:nn)=ppval(ppspline,betacr(i+(nn-1)/2,1:nn));
end
lam_V_qII=zeros(nn,(nn+1)/2);
the_V_qII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
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

M=funfVsII1;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_V_qII*180/pi,the_V_qII*180/pi,Mneg,v,'k--');
contour(lam_V_qII*180/pi,the_V_qII*180/pi,Mpos,v,'k-');
%contour(lam_V_qII*180/pi,the_V_qII*180/pi,funfVsII1,v);
hold on;axis([lmin lmax temin temax umin umax]);

%%     QUADRANT IV DE FACE V

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

M=funfVsIV1;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_V_qIV1*180/pi,the_V_qIV*180/pi,Mneg,v,'k--');
contour(lam_V_qIV1*180/pi,the_V_qIV*180/pi,Mpos,v,'k-');
%contour(lam_V_qIV1*180/pi,the_V_qIV*180/pi,funfVsIV1,v);
hold on;axis([lmin lmax temin temax umin umax]);

%% ************************************************************************
% face VI
%% ************************************************************************

%% QUADRANT I OF FACE VI
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
 lam_VI_a=lam_VI_qI((nn+1)/2:nn,:);
 the_VI_a=the_VI_qI((nn+1)/2:nn,:);
 funfVIa=funfVIsI((nn+1)/2:-1:1,:);

M=funfVIa;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_VI_a*180/pi,the_VI_a*180/pi,Mneg,v,'k--');
contour(lam_VI_a*180/pi,the_VI_a*180/pi,Mpos,v,'k-');
% contour(lam_VI_a*180/pi,the_VI_a*180/pi,funfVIa,v);
 hold on;axis([lmin lmax temin temax umin umax]);

 lam_VI_b=lam_VI_qI(1:(nn+1)/2,:)+2*pi;
 the_VI_b=the_VI_qI(1:(nn+1)/2,:);
 funfVIb=funfVIsI(nn:-1:(nn+1)/2,:);
 
  M=funfVIb;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_VI_b*180/pi,the_VI_b*180/pi,Mneg,v,'k--');
contour(lam_VI_b*180/pi,the_VI_b*180/pi,Mpos,v,'k-');
% contour(lam_VI_b*180/pi,the_VI_b*180/pi,funfVIb,v);
 hold on;axis([lmin lmax temin temax umin umax]);

%% QUADRANT III OF FACE VI 

funfVIsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVIe(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funfVIsIII(1:nn,j)=ppval(ppspline,alfa1(nn:-1:1,j));
end
lam_VI_qIII=zeros(nn,(nn+1)/2);
the_VI_qIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2,
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

 M=funfVIsIII;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour((lam_VI_qIII+pi)*180/pi,the_VI_qIII*180/pi,Mneg,v,'k--');
contour((lam_VI_qIII+pi)*180/pi,the_VI_qIII*180/pi,Mpos,v,'k-');
%contour((lam_VI_qIII+pi)*180/pi,the_VI_qIII*180/pi,funfVIsIII,v);
hold on;axis([lmin lmax temin temax umin umax]);

%%     QUADRANT II DE FACE VI

funfVIsII=zeros((nn+1)/2,nn); % 
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

M=funfVIsII1;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_VI_qII*180/pi,the_VI_qII*180/pi,Mneg,v,'k--');
contour(lam_VI_qII*180/pi,the_VI_qII*180/pi,Mpos,v,'k-');
%contour(lam_VI_qII*180/pi,the_VI_qII*180/pi,funfVIsII1,v);
hold on;axis([lmin lmax temin temax umin umax]);

%%     QUADRANT IV DE FACE VI

funfVIsIV=zeros((nn+1)/2,nn);
for i=1:(nn+1)/2, 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVIe(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVIsIV(i,1:nn)=ppval(ppspline,betacr(i,nn:-1:1));
end
lam_VI_qIV=zeros(nn,(nn+1)/2);
the_VI_qIV=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, 
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
funfVIsIV1=funfVIsIV((nn+1)/2:-1:1,1:nn)';
lam_VI_qIV1=lam_VI_qIV+(2*pi)*(lam_VI_qIV<0);

M=funfVIsIV1;
Mneg=M.*(M<dat);
Mpos=M-Mneg;
contour(lam_VI_qIV1*180/pi,the_VI_qIV*180/pi,Mneg,v,'k--');
contour(lam_VI_qIV1*180/pi,the_VI_qIV*180/pi,Mpos,v,'k-');
%contour(lam_VI_qIV1*180/pi,the_VI_qIV*180/pi,funfVIsIV1,v);
hold on;axis([lmin lmax temin temax umin umax]);


%% ************************************************************************
% TRACE DES FRONTIERES DES PATCHS DE LA CUBED SPHERE GRID
%% ************************************************************************

plot3(lam_Ia((nn+1)/2,:)*180/pi,the_Ia((nn+1)/2,:),funfIea(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ia(:,1)*180/pi,the_Ia(:,1),funfIea(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ia(:,nn)*180/pi,the_Ia(:,nn),funfIea(:,nn)+0.01,'k','LineWidth',1.25); hold on;

plot3(lam_II(1,:)*180/pi,the_II(1,:),funfIIe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(nn,:)*180/pi,the_II(nn,:),funfIIe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,1)*180/pi,the_II(:,1),funfIIe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,nn)*180/pi,the_II(:,nn),funfIIe(:,nn)+0.01,'k','LineWidth',1.25); hold on;

plot3(lam_III(1,:)*180/pi,the_III(1,:),funfIIIe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(nn,:)*180/pi,the_III(nn,:),funfIIIe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(:,1)*180/pi,the_III(:,1),funfIIIe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(:,nn)*180/pi,the_III(:,nn),funfIIIe(:,nn)+0.01,'k','LineWidth',1.25); hold on;

plot3(lam_IV(1,:)*180/pi,the_IV(1,:),funfIVe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(nn,:)*180/pi,the_IV(nn,:),funfIVe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,1)*180/pi,the_IV(:,1),funfIVe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,nn)*180/pi,the_IV(:,nn),funfIVe(:,nn)+0.01,'k','LineWidth',1.25); hold on;

plot3(lam_Ib(1,:)*180/pi,the_Ib(1,:),funfIeb(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ib(:,1)*180/pi,the_Ib(:,1),funfIeb(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ib(:,nn)*180/pi,the_Ib(:,nn),funfIeb(:,nn)+0.01,'k','LineWidth',1.25); hold on;

xlabel('longitude')
ylabel('latitude')

shading interp;
set(gca, 'CLim', [umin, umax]);
view(2);
