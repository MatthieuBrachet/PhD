function []=plot_cs7(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius;
global xi eta dxi deta xx yy delta deltab dga;
global alfa beta;
global alfacr betacr;
global alfa1;
global alfag betag;
 %
 % DESSIN DE TYPE CYLINDRIQUE 
 [lam_I,the_I,rwk]=cart2sph(x_fI,y_fI,z_fI);
 [lam_II,the_II,rwk]=cart2sph(x_fII,y_fII,z_fII);
 [lam_III,the_III,rwk]=cart2sph(x_fIII,y_fIII,z_fIII);
 [lam_IV,the_IV,rwk]=cart2sph(x_fIV,y_fIV,z_fIV);
 [lam_V,the_V,rwk]=cart2sph(x_fV,y_fV,z_fV);
 [lam_VI,the_VI,rwk]=cart2sph(x_fVI,y_fVI,z_fVI);
 %lam_I(:,1)
 %
 lam_I=lam_I+2*pi*(lam_I<-eps); %% TEST A ZERO RECTIFIE.
 lam_II=lam_II+2*pi*(lam_II<-eps);
 lam_III=lam_III+2*pi*(lam_III<-eps);
 lam_IV=lam_IV+2*pi*(lam_IV<-eps);
 lam_V=lam_V+2*pi*(lam_V<-eps);
 lam_VI=lam_VI+2*pi*(lam_VI<-eps);
 %lam_I(:,1)
 %%
lam_Ia=lam_I((nn+1)/2:nn,1:nn);
%max(max(lam_Ia))
%the_Ia=zeros((nn+1)/2,nn);
the_Ia=the_I((nn+1)/2:nn,1:nn);
%
lam_Ib=lam_I(1:(nn+1)/2,1:nn);
lam_Ib((nn+1)/2,:)=lam_Ib((nn+1)/2,:)+2*pi;
the_Ib=the_I(1:(nn+1)/2,1:nn);
%
funfIea=funfIe((nn+1)/2:nn,1:nn);
funfIeb=funfIe(1:(nn+1)/2,1:nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis manual;
lmin=0;lmax=2*pi;temin=-pi/2;temax=pi/2;
% umin=0; umax=1000;
umin=-3; umax=3;
axis([lmin lmax temin temax umin umax]);
%
surf(lam_Ia,the_Ia,funfIea);hold on;axis([lmin lmax temin temax umin umax]);
surf(lam_II,the_II,funfIIe);hold on;axis([lmin lmax temin temax umin umax]);
surf(lam_III,the_III,funfIIIe);hold on;axis([lmin lmax temin temax umin umax]); 
surf(lam_IV,the_IV,funfIVe);hold on;axis([lmin lmax temin temax umin umax]);
surf(lam_Ib,the_Ib,funfIeb);hold on;axis([lmin lmax temin temax umin umax]);
%
% face V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% QUADRANT I OF FACE V
funfVsI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2 % LOOP ON THE eta OF FACE V
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVe(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funfVsI(1:nn,j)=ppval(ppspline,alfa1(1:nn,j));
end
lam_V_qI=zeros(nn,(nn+1)/2);
the_V_qI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %     QUADRANT I DE FACE V
    etawk=eta(j);
    for i=1:nn % FORMULES FACE V
         alfawk=alfa1(i,j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qI(i,j),the_V_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qI(1:nn,(nn+1)/2)=lam_V_qI(1:nn,1);% RECTIFICATION LONGITUDE POUR LE POLE NORD
% 2 TRACES POUR L'ORIGINE DES LAMBDA
% patch a gauche de lambda=0 a pi/4
lam_V_a=lam_V_qI((nn+1)/2:nn,:);
the_V_a=the_V_qI((nn+1)/2:nn,:);
funfVa=funfVsI((nn+1)/2:nn,:);
surf(lam_V_a,the_V_a,funfVa);hold on;axis([lmin lmax temin temax umin umax]);
% patche a droite du dessin de (2*pi)-pi/4 a 2*pi
lam_V_b=lam_V_qI(1:(nn+1)/2,:)+2*pi;
the_V_b=the_V_qI(1:(nn+1)/2,:);
funfVb=funfVsI(1:(nn+1)/2,:);
surf(lam_V_b,the_V_b,funfVb);hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% QUADRANT III OF FACE V ; REVOIR ADRESSAGE: PAS CLAIR
funfVsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, % LOOP ON THE eta OF FACE V
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVe(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    %funfVsIII(1:nn,j)=ppval(ppspline,alfa1(1:nn,nn+1-j)); % JPC
    funfVsIII(1:nn,j)=ppval(ppspline,alfa1(nn:-1:1,nn+1-j));
end
lam_V_qIII=zeros(nn,(nn+1)/2);
the_V_qIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %     QUADRANT III DE FACE V
    etawk=eta(j);
    for i=1:nn % FORMULES FACE V
         alfawk=alfa1(i,nn+1-j);
         xx_wk1 = -cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = cos(alfawk)*cos(etawk);
         [lam_V_qIII(i,j),the_V_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIII(1:nn,(nn+1)/2)=lam_V_qIII(1:nn,1);% RECTIFICATION LONGITUDE POUR LE POLE NORD
surf(lam_V_qIII+pi,the_V_qIII,funfVsIII);
%colorbar;
hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT II DE FACE V
funfVsII=zeros((nn+1)/2,nn); % 
for i=(nn+1)/2:-1:1, % LOOP ON THE xi OF FACE V
    betaspline(1:nn)=beta(i+(nn-1)/2,1:nn);
    funspline(1:nn)=funfVe(i+(nn-1)/2,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVsII(i,1:nn)=ppval(ppspline,betacr(i+(nn-1)/2,1:nn));
end
lam_V_qII=zeros(nn,(nn+1)/2);
the_V_qII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %     QUADRANT II DE FACE V
    xiwk=xi(nn-j+1);
    for i=1:nn % FORMULES FACE V
         betawk = -betacr(j,i); % CONTROLLER
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qII(i,j),the_V_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qII(1:nn,(nn+1)/2)=lam_V_qII(1:nn,1);% RECTIFICATION LONGITUDE POUR LE POLE NORD
funfVsII1=funfVsII((nn+1)/2:-1:1,:)';
surf(lam_V_qII,the_V_qII,funfVsII1);hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT IV DE FACE V
funfVsIV=zeros((nn+1)/2,nn); % 
for i=1:(nn+1)/2, % LOOP ON THE xi OF FACE V
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVe(i,1:nn);
    ppspline=spline(betaspline,funspline);
%     funfVsIV(i,1:nn)=ppval(ppspline,betacr(i,1:nn)); % JPC
     funfVsIV(i,1:nn)=ppval(ppspline,betacr(i,nn:-1:1));
end
lam_V_qIV=zeros(nn,(nn+1)/2);
the_V_qIV=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %  QUADRANT IV DE FACE V
    xiwk=xi(j);
    for i=1:nn % FORMULES FACE V
         betawk = betacr(j,nn-i+1); % CONTROLLER
         xx_wk1 = -sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = cos(betawk)*cos(-xiwk);
         [lam_V_qIV(i,j),the_V_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_V_qIV(1:nn,(nn+1)/2)=lam_V_qIV(1:nn,1); % RECTIFICATION LONGITUDE POUR LE POLE NORD
funfVsIV1=funfVsIV(1:(nn+1)/2,:)';
lam_V_qIV1=lam_V_qIV+(2*pi)*(lam_V_qIV<0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf(lam_V_qIV1,the_V_qIV,funfVsIV1);hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIN TRACE FACE V  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE VI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% QUADRANT I OF FACE VI
funfVIsI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2 % LOOP ON THE eta OF FACE V
    alfaspline(1:nn)=alfa(1:nn,j+(nn+1)/2-1);
    funspline(1:nn)=funfVIe(1:nn,j+(nn+1)/2-1);
    ppspline=spline(alfaspline,funspline);
    funfVIsI(1:nn,j)=ppval(ppspline,alfa1(1:nn,j+(nn+1)/2-1));
end
lam_VI_qI=zeros(nn,(nn+1)/2);
the_VI_qI=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, 
%     etawk=eta(nn-j+1); % JPC
    etawk=eta(j+(nn+1)/2-1); 
    for i=1:nn % FORMULES FACE VI
         alfawk=alfa1(i,(nn+1)/2-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qI(i,j),the_VI_qI(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
% lam_VI_qI(1:nn,(nn+1)/2)=lam_VI_qI(1:nn,1);
% LE POLE NORD % JPC
lam_VI_qI(1:nn,1)=lam_VI_qI(1:nn,(nn+1)/2);
%max(lam_VI_qI)
%surf(lam_VI_qI,the_VI_qI,funfVIsI(nn:-1:1,:));hold on;
% % 2 TRACES POUR L'ORIGINE DES LAMBDA
% % patch a gauche de lambda=0 a pi/4
 lam_VI_a=lam_VI_qI((nn+1)/2:nn,:);
 the_VI_a=the_VI_qI((nn+1)/2:nn,:);
%  funfVIa=funfVIsI((nn+1)/2:nn,(nn+1)/2:-1:1); % JPC
 funfVIa=funfVIsI((nn+1)/2:-1:1,:);
 %funfVIa=funfVIsI(nn:-1:(nn+1)/2,:);
 surf(lam_VI_a,the_VI_a,funfVIa);hold on;axis([lmin lmax temin temax umin umax]);
% % patche a droite du dessin de (2*pi)-pi/4 a 2*pi
 lam_VI_b=lam_VI_qI(1:(nn+1)/2,:)+2*pi;
 the_VI_b=the_VI_qI(1:(nn+1)/2,:);
 funfVIb=funfVIsI(nn:-1:(nn+1)/2,:);
 %funfVIb=funfVIsI((nn+1)/2:-1:1,:);
 surf(lam_VI_b,the_VI_b,funfVIb);hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% QUADRANT III OF FACE VI ; REVOIR ADRESSAGE: PAS CLAIR
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
surf(lam_VI_qIII+pi,the_VI_qIII,funfVIsIII);
hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT II DE FACE VI
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
surf(lam_VI_qII,the_VI_qII,funfVIsII1);hold on;axis([lmin lmax temin temax umin umax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT IV DE FACE VI
funfVIsIV=zeros((nn+1)/2,nn); % 
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
         betawk = betacr((nn+1)/2-j+1,nn-i+1); % CONTROLLER
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qIV(i,j),the_VI_qIV(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIV(1:nn,1)=lam_VI_qIV(1:nn,(nn+1)/2); % RECTIFICATION LONGITUDE POUR LE POLE NORD
funfVIsIV1=funfVIsIV((nn+1)/2:-1:1,1:nn)';
lam_VI_qIV1=lam_VI_qIV+(2*pi)*(lam_VI_qIV<0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf(lam_VI_qIV1,the_VI_qIV,funfVIsIV1);hold on;axis([lmin lmax temin temax umin umax]);
view(3);
%% ************************************************************************
% TRACE DES FRONTIERES DES PATCHS DE LA CUBED SPHERE GRID
% DESSIN DE TYPE CYLINDRIQUE ====
plot3(lam_Ia((nn+1)/2,:),the_Ia((nn+1)/2,:),funfIea(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ia(:,1),the_Ia(:,1),funfIea(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ia(:,nn),the_Ia(:,nn),funfIea(:,nn)+0.01,'k','LineWidth',1.25); hold on;
%
plot3(lam_II(1,:),the_II(1,:),funfIIe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(nn,:),the_II(nn,:),funfIIe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,1),the_II(:,1),funfIIe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_II(:,nn),the_II(:,nn),funfIIe(:,nn)+0.01,'k','LineWidth',1.25); hold on;
%
plot3(lam_III(1,:),the_III(1,:),funfIIIe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(nn,:),the_III(nn,:),funfIIIe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(:,1),the_III(:,1),funfIIIe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_III(:,nn),the_III(:,nn),funfIIIe(:,nn)+0.01,'k','LineWidth',1.25); hold on;
%
plot3(lam_IV(1,:),the_IV(1,:),funfIVe(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(nn,:),the_IV(nn,:),funfIVe(nn,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,1),the_IV(:,1),funfIVe(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_IV(:,nn),the_IV(:,nn),funfIVe(:,nn)+0.01,'k','LineWidth',1.25); hold on;
%
plot3(lam_Ib(1,:),the_Ib(1,:),funfIeb(1,:)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ib(:,1),the_Ib(:,1),funfIeb(:,1)+0.01,'k','LineWidth',1.25); hold on;
plot3(lam_Ib(:,nn),the_Ib(:,nn),funfIeb(:,nn)+0.01,'k','LineWidth',1.25); hold on;
