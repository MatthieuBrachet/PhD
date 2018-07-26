function []=plot_cs9(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
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

% clear all; 
% % load ('exec68/NL08-Jan-09-2013a');
% % load ('exec68/NL16-Jan-09-2013a');
% % load ('exec68/NL32-Jan-09-2013a');
%  load ('exec68/NL64-Jan-09-2013a');
% %  load ('exec58/March-26-2012b');
% %  load ('exec58/March-26-2012c');
% % load ('exec58/March-26-2012d');
%  figure(1);
%  plot_cs5(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe);colorbar;
%  %%%%%%
%  figure(2);
%  plot_cs5(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI);colorbar;
%  %%%%%%%%
%  figure(3);
%   plot_cs5(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI);colorbar;
%   %%%%%%%%%%%%
%  figure(4);
%  plot(xsecs,er1,'k-');hold on;grid;
%  plot(xsecs,er2,'k--');hold on;
%  plot(xsecs,erinfty,'k.');
%  %
%  figure(5);
%  plot(xsecs,ermax,'k-');hold on;grid;
%  plot(xsecs,ermin,'k--');
 %
 % DESSIN DE TYPE CYLINDRIQUE 
 [lam_I,the_I,rwk]=cart2sph(x_fI,y_fI,z_fI);
 [lam_II,the_II,rwk]=cart2sph(x_fII,y_fII,z_fII);
 [lam_III,the_III,rwk]=cart2sph(x_fIII,y_fIII,z_fIII);
 [lam_IV,the_IV,rwk]=cart2sph(x_fIV,y_fIV,z_fIV);
 [lam_V,the_V,rwk]=cart2sph(x_fV,y_fV,z_fV);
 [lam_VI,the_VI,rwk]=cart2sph(x_fVI,y_fVI,z_fVI);
 %
 lam_I=lam_I+2*pi*(lam_I<-eps);
 lam_II=lam_II+2*pi*(lam_II<-eps);
 lam_III=lam_III+2*pi*(lam_III<-eps);
 lam_IV=lam_IV+2*pi*(lam_IV<-eps);
 lam_V=lam_V+2*pi*(lam_V<-eps);
 lam_VI=lam_VI+2*pi*(lam_VI<-eps);
 %%
lam_Ia=lam_I((nn+1)/2:nn,1:nn);
the_Ia=zeros((nn+1)/2,1:nn);
the_Ia=the_I((nn+1)/2:nn,1:nn);
%
lam_Ib=lam_I(1:(nn+1)/2,1:nn);
lam_Ib((nn+1)/2,:)=lam_Ib((nn+1)/2,:)+2*pi;
the_Ib=the_I(1:(nn+1)/2,1:nn);
%
funfIea=funfIe((nn+1)/2:nn,1:nn);
funfIeb=funfIe(1:(nn+1)/2,1:nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contour(lam_Ia,the_Ia,funfIea); hold on;
% contour(lam_II,the_II,funfIIe); hold on;
% contour(lam_III,the_III,funfIIIe);hold on; 
% contour(lam_IV,the_IV,funfIVe);hold on; 
% contour(lam_Ib,the_Ib,funfIeb); hold on;
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
%min(min(lam_V_qI))
% 2 TRACES POUR L'ORIGINE DES LAMBDA
% patch a gauche de lambda=0 a pi/4
lam_V_a=lam_V_qI((nn+1)/2:nn,:);
the_V_a=the_V_qI((nn+1)/2:nn,:);
funfVa=funfVsI((nn+1)/2:nn,:);
% contour(lam_V_a,the_V_a,funfVa);hold on;
% patche a droite du dessin de (2*pi)-pi/4 a 2*pi
lam_V_b=lam_V_qI(1:(nn+1)/2,:)+2*pi;
the_V_b=the_V_qI(1:(nn+1)/2,:);
funfVb=funfVsI(1:(nn+1)/2,:);
% contour(lam_V_b,the_V_b,funfVb);hold on;
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
% contour(lam_V_qIII+pi,the_V_qIII,funfVsIII);
colorbar;hold on;
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
% contour(lam_V_qII,the_V_qII,funfVsII1);hold on;
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
% contour(lam_V_qIV1,the_V_qIV,funfVsIV1);hold on;
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
for j=1:(nn+1)/2, %     QUADRANT I DE FACE VI
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
% lam_VI_qI(1:nn,(nn+1)/2)=lam_VI_qI(1:nn,1);% RECTIFICATION LONGITUDE POUR
% LE POLE NORD % JPC
lam_VI_qI(1:nn,1)=lam_VI_qI(1:nn,(nn+1)/2);% RECTIFICATION LONGITUDE POUR LE POLE NORD
%contour(lam_VI_qI,the_VI_qI,funfVIsI(nn:-1:1,:));hold on;
% % 2 TRACES POUR L'ORIGINE DES LAMBDA
% % patch a gauche de lambda=0 a pi/4
 lam_VI_a=lam_VI_qI((nn+1)/2:nn,:);
 the_VI_a=the_VI_qI((nn+1)/2:nn,:);
%  funfVIa=funfVIsI((nn+1)/2:nn,(nn+1)/2:-1:1); % JPC
 funfVIa=funfVIsI((nn+1)/2:-1:1,:);
 %funfVIa=funfVIsI(nn:-1:(nn+1)/2,:);
%  contour(lam_VI_a,the_VI_a,funfVIa);hold on;
% % patche a droite du dessin de (2*pi)-pi/4 a 2*pi
 lam_VI_b=lam_VI_qI(1:(nn+1)/2,:)+2*pi;
 the_VI_b=the_VI_qI(1:(nn+1)/2,:);
 funfVIb=funfVIsI(nn:-1:(nn+1)/2,:);
 %funfVIb=funfVIsI((nn+1)/2:-1:1,:);
%  contour(lam_VI_b,the_VI_b,funfVIb);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% QUADRANT III OF FACE VI ; REVOIR ADRESSAGE: PAS CLAIR
funfVIsIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, % LOOP ON THE eta OF FACE VI
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVIe(1:nn,j);
    ppspline=spline(alfaspline,funspline);
  %  funfVIsIII(1:nn,j)=ppval(ppspline,alfa1(1:nn,j)); % JPC
    funfVIsIII(1:nn,j)=ppval(ppspline,alfa1(nn:-1:1,j));
end
lam_VI_qIII=zeros(nn,(nn+1)/2);
the_VI_qIII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %     QUADRANT III DE FACE VI
    etawk=eta(nn-j+1);
    for i=1:nn % FORMULES FACE VI
         alfawk=alfa1(nn-i+1,nn-j+1);
         xx_wk1 = cos(alfawk)*sin(etawk);
         yy_wk1 = sin(alfawk);   
         zz_wk1 = -cos(alfawk)*cos(etawk);
         [lam_VI_qIII(i,j),the_VI_qIII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qIII(1:nn,(nn+1)/2)=lam_VI_qIII(1:nn,1);% RECTIFICATION LONGITUDE POUR LE POLE SUD
% contour(lam_VI_qIII+pi,the_VI_qIII,funfVIsIII);
colorbar;hold on;
% break;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT II DE FACE VI
funfVIsII=zeros((nn+1)/2,nn); % 
for i=(nn+1)/2:-1:1, % LOOP ON THE xi OF FACE VI
    betaspline(1:nn)=beta(i+(nn-1)/2,1:nn);
    funspline(1:nn)=funfVIe(i+(nn-1)/2,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVIsII(i,1:nn)=ppval(ppspline,betacr(i+(nn-1)/2,1:nn));
end
lam_VI_qII=zeros(nn,(nn+1)/2);
the_VI_qII=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %     QUADRANT II DE FACE VI
    xiwk=xi(nn-j+1);
    for i=1:nn % FORMULES FACE VI
         betawk = betacr(j,i); % CONTROLLER
         xx_wk1 = sin(betawk);
         yy_wk1 = -cos(betawk)*sin(-xiwk);   
         zz_wk1 = -cos(betawk)*cos(-xiwk);
         [lam_VI_qII(i,j),the_VI_qII(i,j),rrwk] = cart2sph(xx_wk1,yy_wk1,zz_wk1);
    end
end
lam_VI_qII(1:nn,(nn+1)/2)=lam_VI_qII(1:nn,1);% RECTIFICATION LONGITUDE POUR LE POLE NORD
funfVIsII1=funfVIsII((nn+1)/2:-1:1,nn:-1:1)';
% contour(lam_VI_qII,the_VI_qII,funfVIsII1);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     QUADRANT IV DE FACE VI
funfVIsIV=zeros((nn+1)/2,nn); % 
for i=1:(nn+1)/2, % LOOP ON THE xi OF FACE VI
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVIe(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funfVIsIV(i,1:nn)=ppval(ppspline,betacr(i,nn:-1:1));
end
lam_VI_qIV=zeros(nn,(nn+1)/2);
the_VI_qIV=zeros(nn,(nn+1)/2);
for j=1:(nn+1)/2, %  QUADRANT IV DE FACE VI
    xiwk=xi((nn+1)/2-j+1);
    for i=1:nn % FORMULES FACE VI
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
% contour(lam_VI_qIV1,the_VI_qIV,funfVIsIV1);hold on;
shading interp;
axis([0 2*pi -pi/2 pi/2]);
% figure(13);
lam_g=[lam_VI_a,lam_Ia,lam_V_a;
    lam_VI_qII,lam_II,lam_V_qII;
    lam_VI_qIII+pi,lam_III,lam_V_qIII+pi;
    lam_VI_qIV1,lam_IV,lam_V_qIV1;
    lam_VI_b,lam_Ib,lam_V_b];
the_g=[the_VI_a,the_Ia,the_V_a;
    the_VI_qII,the_II,the_V_qII;
    the_VI_qIII,the_III,the_V_qIII;
    the_VI_qIV,the_IV,the_V_qIV;
    the_VI_b,the_Ib,the_V_b];
fun_g=[funfVIa,funfIea,funfVa;
       funfVIsII1,funfIIe,funfVsII1;
       funfVIsIII,funfIIIe,funfVsIII;
       funfVIsIV1,funfIVe,funfVsIV1;
       funfVIb,funfIeb,funfVb];
   size(lam_g)
   size(the_g)
   size(fun_g)
contour (lam_g,the_g,fun_g);
% contour(lam_Ia,the_Ia,funfIea); hold on;
% contour(lam_II,the_II,funfIIe); hold on;
% contour(lam_III,the_III,funfIIIe);hold on; 
% contour(lam_IV,the_IV,funfIVe);hold on; 
% contour(lam_Ib,the_Ib,funfIeb); hold on;
% 
% contour(lam_V_a,the_V_a,funfVa);hold on;
% contour(lam_V_qII,the_V_qII,funfVsII1);hold on;
% contour(lam_V_qIII+pi,the_V_qIII,funfVsIII);
% contour(lam_V_qIV1,the_V_qIV,funfVsIV1);hold on;
% contour(lam_V_b,the_V_b,funfVb);hold on;
% 
% contour(lam_VI_a,the_VI_a,funfVIa);hold on;
% contour(lam_VI_qII,the_VI_qII,funfVIsII1);hold on;
% contour(lam_VI_qIII+pi,the_VI_qIII,funfVIsIII);
% contour(lam_VI_qIV1,the_VI_qIV,funfVIsIV1);hold on;
% contour(lam_VI_b,the_VI_b,funfVIb);hold on;

% lam_III(1:nn,1)'
% (lam_V_qI(1:nn,1)+pi)'
% funfIIIe(1:nn,(nn+1)/2)'
% funfVsIII(1:nn,1)'
% lam_V_====
% funfIIIe(1:nn,1)'
% funfVIsIII(1:nn,(nn+1)/2)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIN TRACE FACE VI  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % axis square;