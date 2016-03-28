% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
clear all;clc; close all;
global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius dxi
global k1 omega;
global ite
global u0;
global t_final;

%% space data
n=31;
mod72; 

%% time data
tinit=0;
tmax=t_final;
crk=2.5;
stab=0.25;
ddt=stab*crk/(2*omega);
ntime=2;
itemax=ntime;
cfl=(ddt*u0)/(radius*dxi);

%% initial condition
vit1=fun8(x_fI,y_fI,z_fI,0);


[mfunfI]=fun8(x_fI,y_fI,z_fI,tinit); 
[mfunfII]=fun8(x_fII,y_fII,z_fII,tinit);
[mfunfIII]=fun8(x_fIII,y_fIII,z_fIII,tinit);
[mfunfIV]=fun8(x_fIV,y_fIV,z_fIV,tinit);
[mfunfV]=fun8(x_fV,y_fV,z_fV,tinit);
[mfunfVI]=fun8(x_fVI,y_fVI,z_fVI,tinit);

mfunfInew=zeros(nn,nn,4);mfunfIInew=zeros(nn,nn,4);mfunfIIInew=zeros(nn,nn,4);
mfunfIVnew=zeros(nn,nn,4);mfunfVnew=zeros(nn,nn,4);mfunfVInew=zeros(nn,nn,4);
kk0_I=zeros(nn,nn,4);kk1_I=zeros(nn,nn,4);kk2_I=zeros(nn,nn,4);kk3_I=zeros(nn,nn,4);
kk0_II=zeros(nn,nn,4);kk1_II=zeros(nn,nn,4);kk2_II=zeros(nn,nn,4);kk3_II=zeros(nn,nn,4);
kk0_III=zeros(nn,nn,4);kk1_III=zeros(nn,nn,4);kk2_III=zeros(nn,nn,4);kk3_III=zeros(nn,nn,4);
kk0_IV=zeros(nn,nn,4);kk1_IV=zeros(nn,nn,4);kk2_IV=zeros(nn,nn,4);kk3_IV=zeros(nn,nn,4);
kk0_V=zeros(nn,nn,4);kk1_V=zeros(nn,nn,4);kk2_V=zeros(nn,nn,4);kk3_V=zeros(nn,nn,4);
kk0_VI=zeros(nn,nn,4);kk1_VI=zeros(nn,nn,4);kk2_VI=zeros(nn,nn,4);kk3_VI=zeros(nn,nn,4);
mfunftI=zeros(nn,nn,4); mfunftII=zeros(nn,nn,4); mfunftIII=zeros(nn,nn,4);
mfunftIV=zeros(nn,nn,4); mfunftV=zeros(nn,nn,4); mfunftVI=zeros(nn,nn,4);


time=tinit;
er1=zeros(1,itemax);er2=zeros(1,itemax);erinfty=zeros(1,itemax);ermax=zeros(1,itemax);ermin=zeros(1,itemax);
xdays=zeros(1,itemax);
integral=zeros(1,itemax);
for ite=1:itemax
ite
%% *** filtrage ***********************************************************
for k1=1:4,
[mfunftI(:,:,k1),mfunftII(:,:,k1),mfunftIII(:,:,k1),mfunftIV(:,:,k1),mfunftV(:,:,k1),mfunftVI(:,:,k1)]=...
    ftr72(mfunfI(:,:,k1),mfunfII(:,:,k1),mfunfIII(:,:,k1),mfunfIV(:,:,k1),mfunfV(:,:,k1),mfunfVI(:,:,k1),n,nn);
end
mfunfI=mfunftI;mfunfII=mfunftII;mfunfIII=mfunftIII;
mfunfIV=mfunftIV;mfunfV=mfunftV;mfunfVI=mfunftVI;

%% *** Iterations RK4 *****************************************************
%% CALCUL KK0
[mflu_I,mflu_II,mflu_III,mflu_IV,mflu_V,mflu_VI]=...
    flu72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

%[mfor_I,mfor_II,mfor_III,mfor_IV,mfor_V,mfor_VI]=...
%     for72(mfun2_fI,mfun2_fII,mfun2_fIII,mfun2_fIV,mfun2_fV,mfun2_fVI,time,n,nn);

kk0_I(1:nn,1:nn,1:4) = -mflu_I(1:nn,1:nn,1:4); %+  mfor_I(1:nn,1:nn,1:4);
kk0_II(1:nn,1:nn,1:4) = -mflu_II(1:nn,1:nn,1:4);% +  mfor_II(1:nn,1:nn,1:4);
kk0_III(1:nn,1:nn,1:4) = -mflu_III(1:nn,1:nn,1:4);% +  mfor_III(1:nn,1:nn,1:4);
kk0_IV(1:nn,1:nn,1:4) = -mflu_IV(1:nn,1:nn,1:4); %+  mfor_IV(1:nn,1:nn,1:4);
kk0_V(1:nn,1:nn,1:4) = -mflu_V(1:nn,1:nn,1:4); %+  mfor_V(1:nn,1:nn,1:4);
kk0_VI(1:nn,1:nn,1:4) = -mflu_VI(1:nn,1:nn,1:4); %+  mfor_VI(1:nn,1:nn,1:4);

%% CALCUL KK1
mfun1_fI=mfunfI+0.5*ddt*kk0_I;
mfun1_fII=mfunfII+0.5*ddt*kk0_II;
mfun1_fIII=mfunfIII+0.5*ddt*kk0_III;
mfun1_fIV=mfunfIV+0.5*ddt*kk0_IV;
mfun1_fV=mfunfV+0.5*ddt*kk0_V;
mfun1_fVI=mfunfVI+0.5*ddt*kk0_VI;

[mflu_I,mflu_II,mflu_III,mflu_IV,mflu_V,mflu_VI]=...
    flu72(mfun1_fI,mfun1_fII,mfun1_fIII,mfun1_fIV,mfun1_fV,mfun1_fVI,n,nn);

kk1_I(1:nn,1:nn,1:4) = -mflu_I(1:nn,1:nn,1:4); % +  mfor_I(1:nn,1:nn,1:4);
kk1_II(1:nn,1:nn,1:4) = -mflu_II(1:nn,1:nn,1:4); % +  mfor_II(1:nn,1:nn,1:4);
kk1_III(1:nn,1:nn,1:4) = -mflu_III(1:nn,1:nn,1:4);% +  mfor_III(1:nn,1:nn,1:4);
kk1_IV(1:nn,1:nn,1:4) = -mflu_IV(1:nn,1:nn,1:4); %+  mfor_IV(1:nn,1:nn,1:4);
kk1_V(1:nn,1:nn,1:4) = -mflu_V(1:nn,1:nn,1:4); %+  mfor_V(1:nn,1:nn,1:4);
kk1_VI(1:nn,1:nn,1:4) = -mflu_VI(1:nn,1:nn,1:4); %+  mfor_VI(1:nn,1:nn,1:4);

%% CALCUL KK2
mfun2_fI=mfunfI+0.5*ddt*kk1_I;
mfun2_fII=mfunfII+0.5*ddt*kk1_II;
mfun2_fIII=mfunfIII+0.5*ddt*kk1_III;
mfun2_fIV=mfunfIV+0.5*ddt*kk1_IV;
mfun2_fV=mfunfV+0.5*ddt*kk1_V;
mfun2_fVI=mfunfVI+0.5*ddt*kk1_VI;

[mflu_I,mflu_II,mflu_III,mflu_IV,mflu_V,mflu_VI]=...
    flu72(mfun2_fI,mfun2_fII,mfun2_fIII,mfun2_fIV,mfun2_fV,mfun2_fVI,n,nn);

kk2_I(1:nn,1:nn,1:4) = -mflu_I(1:nn,1:nn,1:4);% +  mfor_I(1:nn,1:nn,1:4);
kk2_II(1:nn,1:nn,1:4) = -mflu_II(1:nn,1:nn,1:4);% +  mfor_II(1:nn,1:nn,1:4);
kk2_III(1:nn,1:nn,1:4) = -mflu_III(1:nn,1:nn,1:4); %+  mfor_III(1:nn,1:nn,1:4);
kk2_IV(1:nn,1:nn,1:4) = -mflu_IV(1:nn,1:nn,1:4); %+  mfor_IV(1:nn,1:nn,1:4);
kk2_V(1:nn,1:nn,1:4) = -mflu_V(1:nn,1:nn,1:4); %+  mfor_V(1:nn,1:nn,1:4);
kk2_VI(1:nn,1:nn,1:4) = -mflu_VI(1:nn,1:nn,1:4); %+  mfor_VI(1:nn,1:nn,1:4);

%% CALCUL KK3
mfun3_fI=mfunfI+ddt*kk2_I;
mfun3_fII=mfunfII+ddt*kk2_II;
mfun3_fIII=mfunfIII+ddt*kk2_III;
mfun3_fIV=mfunfIV+ddt*kk2_IV;
mfun3_fV=mfunfV+ddt*kk2_V;
mfun3_fVI=mfunfVI+ddt*kk2_VI;

[mflu_I,mflu_II,mflu_III,mflu_IV,mflu_V,mflu_VI]=...
    flu72(mfun3_fI,mfun3_fII,mfun3_fIII,mfun3_fIV,mfun3_fV,mfun3_fVI,n,nn);

kk3_I(1:nn,1:nn,1:4) = -mflu_I(1:nn,1:nn,1:4); %+  mfor_I(1:nn,1:nn,1:4);
kk3_II(1:nn,1:nn,1:4) = -mflu_II(1:nn,1:nn,1:4); %+  mfor_II(1:nn,1:nn,1:4);
kk3_III(1:nn,1:nn,1:4) = -mflu_III(1:nn,1:nn,1:4); %+  mfor_III(1:nn,1:nn,1:4);
kk3_IV(1:nn,1:nn,1:4) = -mflu_IV(1:nn,1:nn,1:4); %+  mfor_IV(1:nn,1:nn,1:4);
kk3_V(1:nn,1:nn,1:4) = -mflu_V(1:nn,1:nn,1:4); %+  mfor_V(1:nn,1:nn,1:4);
kk3_VI(1:nn,1:nn,1:4) = -mflu_VI(1:nn,1:nn,1:4); %+  mfor_VI(1:nn,1:nn,1:4);

%% ASSEMBLAGE RK4
mfunfInew(1:nn,1:nn,1:4)=mfunfI(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_I(1:nn,1:nn,1:4)+(1/3)*kk1_I(1:nn,1:nn,1:4)+(1/3)*kk2_I(1:nn,1:nn,1:4)+(1/6)*kk3_I(1:nn,1:nn,1:4));
mfunfIInew(1:nn,1:nn,1:4)=mfunfII(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_II(1:nn,1:nn,1:4)+(1/3)*kk1_II(1:nn,1:nn,1:4)+(1/3)*kk2_II(1:nn,1:nn,1:4)+(1/6)*kk3_II(1:nn,1:nn,1:4));
mfunfIIInew(1:nn,1:nn,1:4)=mfunfIII(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_III(1:nn,1:nn,1:4)+(1/3)*kk1_III(1:nn,1:nn,1:4)+(1/3)*kk2_III(1:nn,1:nn,1:4)+(1/6)*kk3_III(1:nn,1:nn,1:4));
mfunfIVnew(1:nn,1:nn,1:4)=mfunfIV(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_IV(1:nn,1:nn,1:4)+(1/3)*kk1_IV(1:nn,1:nn,1:4)+(1/3)*kk2_IV(1:nn,1:nn,1:4)+(1/6)*kk3_IV(1:nn,1:nn,1:4));
mfunfVnew(1:nn,1:nn,1:4)=mfunfV(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_V(1:nn,1:nn,1:4)+(1/3)*kk1_V(1:nn,1:nn,1:4)+(1/3)*kk2_V(1:nn,1:nn,1:4)+(1/6)*kk3_V(1:nn,1:nn,1:4));
mfunfVInew(1:nn,1:nn,1:4)=mfunfVI(1:nn,1:nn,1:4)+ddt*...
    ((1/6)*kk0_VI(1:nn,1:nn,1:4)+(1/3)*kk1_VI(1:nn,1:nn,1:4)+(1/3)*kk2_VI(1:nn,1:nn,1:4)+(1/6)*kk3_VI(1:nn,1:nn,1:4));

%% *** demi sommes ********************************************************
[mfunfInew,mfunfIInew,mfunfIIInew,mfunfIVnew,mfunfVnew,mfunfVInew]=...
    ds72(mfunfInew,mfunfIInew,mfunfIIInew,mfunfIVnew,mfunfVnew,mfunfVInew,n,nn);

%% *** mise a jour ********************************************************
time=time+ddt;

mfunfI=mfunfInew;mfunfII=mfunfIInew;mfunfIII=mfunfIIInew;
mfunfIV=mfunfIVnew;mfunfV=mfunfVnew;mfunfVI=mfunfVInew;

% %% HISTORIQUE ERREUR
% xdays(ite)=(time-tinit)/(24*3600);
% xsecs(ite)=time-tinit;
% 
% mfunfIe=fun8(x_fI,y_fI,z_fI,time);mfunfIIe=fun8(x_fII,y_fII,z_fII,time);
% mfunfIIIe=fun8(x_fIII,y_fIII,z_fIII,time); mfunfIVe=fun8(x_fIV,y_fIV,z_fIV,time);
% mfunfVe=fun8(x_fV,y_fV,z_fV,time); mfunfVIe=fun8(x_fVI,y_fVI,z_fVI,time);
% err_fI=mfunfI-mfunfIe;err_fII=mfunfII-mfunfIIe;
% err_fIII=mfunfIII-mfunfIIIe;err_fIV=mfunfIV-mfunfIVe;
% err_fV=mfunfV-mfunfVe;err_fVI=mfunfVI-mfunfVIe;
% 
%  str='1';
% [nrmerI(1:4),nrmerII(1:4),nrmerIII(1:4),nrmerIV(1:4),nrmerV(1:4),nrmerVI(1:4),nrmger(1:4)]=...
%     nrm72(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
% [nrmeI(1:4),nrmeII(1:4),nrmeIII(1:4),nrmeIV(1:4),nrmeV(1:4),nrmeVI(1:4),nrmge(1:4)]=...
%     nrm72(mfunfIe,mfunfIIe,mfunfIIIe,mfunfIVe,mfunfVe,mfunfVIe,n,nn,str);
% er1(ite,1:4)=nrmger(1:4)./nrmge(1:4);
% 
% str='2';
% [nrmerI(1:4),nrmerII(1:4),nrmerIII(1:4),nrmerIV(1:4),nrmerV(1:4),nrmerVI(1:4),nrmger(1:4)]=...
%     nrm72(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
% [nrmeI(1:4),nrmeII(1:4),nrmeIII(1:4),nrmeIV(1:4),nrmeV(1:4),nrmeVI(1:4),nrmge(1:4)]=...
%     nrm72(mfunfIe,mfunfIIe,mfunfIIIe,mfunfIVe,mfunfVe,mfunfVIe,n,nn,str);
% er2(ite,1:4)=nrmger./nrmge;
% 
% str='infty';
% [nrmerI(1:4),nrmerII(1:4),nrmerIII(1:4),nrmerIV(1:4),nrmerV(1:4),nrmerVI(1:4),nrmger(1:4)]=...
%     nrm72(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
% [nrmeI(1:4),nrmeII(1:4),nrmeIII(1:4),nrmeIV(1:4),nrmeV(1:4),nrmeVI(1:4),nrmge(1:4)]=...
%     nrm72(mfunfIe,mfunfIIe,mfunfIIIe,mfunfIVe,mfunfVe,mfunfVIe,n,nn,str);
% erinfty(ite,1:4)=nrmger./nrmge;
end

%% *** bilan **************************************************************
disp('1');
er1(itemax)
disp('2');
er2(itemax)
disp('infty');
erinfty(itemax)

%% bilan solution temps final
mfunfIe=fun8(x_fI,y_fI,z_fI,time);mfunfIIe=fun8(x_fII,y_fII,z_fII,time);
mfunfIIIe=fun8(x_fIII,y_fIII,z_fIII,time); mfunfIVe=fun8(x_fIV,y_fIV,z_fIV,time);
mfunfVe=fun8(x_fV,y_fV,z_fV,time); mfunfVIe=fun8(x_fVI,y_fVI,z_fVI,time);
err_fI=mfunfI-mfunfIe;err_fII=mfunfII-mfunfIIe;
err_fIII=mfunfIII-mfunfIIIe;err_fIV=mfunfIV-mfunfIVe;
err_fV=mfunfV-mfunfVe;err_fVI=mfunfVI-mfunfVIe;


% figure(1)
% semilogy(xdays, er1,xdays,er2,xdays,erinfty)
% legend('norm 1','norm 2','norm infty')

errI_38=max(max(abs(err_fI(:,:,4)))); % ERREUR ETA
errII_38=max(max(abs(err_fII(:,:,4)))); % ERREUR ETA
errIII_38=max(max(abs(err_fIII(:,:,4)))); % ERREUR ETA
errIV_38=max(max(abs(err_fIV(:,:,4)))); % ERREUR ETA
errV_38=max(max(abs(err_fV(:,:,4)))); % ERREUR ETA
errVI_38=max(max(abs(err_fVI(:,:,4)))); % ERREUR ETA
err_38=max([errI_38,errII_38,errIII_38,errIV_38,errV_38,errVI_38])

figure(33);
plot_cs5(n,nn,err_fI(:,:,4),err_fII(:,:,4),err_fIII(:,:,4)...
     ,err_fIV(:,:,4),err_fV(:,:,4),err_fVI(:,:,4));colorbar;
title('erreur')

figure(35);
plot_cs5(n,nn,mfunfIe(:,:,4),mfunfIIe(:,:,4),mfunfIIIe(:,:,4)...
     ,mfunfIVe(:,:,4),mfunfVe(:,:,4),mfunfVIe(:,:,4));colorbar;
title('exacte')

figure(37);
 plot_cs5(n,nn,mfunfI(:,:,4),mfunfII(:,:,4),mfunfIII(:,:,4)...
    ,mfunfIV(:,:,4),mfunfV(:,:,4),mfunfVI(:,:,4));colorbar;
title('approchee')