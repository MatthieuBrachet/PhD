%% ************************************************************************
% test for laplacian on the CS.
% *** options :
% scheme : numerical spatial scheme used. 
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or implicit (visbal) filtering.
% alfa_ftr : parameter for implicit fliter (type visbal only, 
%            if alpha_ftr=0, the filter is equivalent to Redonnet filter,
%            if alpha_ftr=0.5, the filter is inexistant).
%
%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr alfa_ftr

opt_ftr='redonnet10';
alfa_ftr=0;
scheme='compact4';

n=31; % for snapshot, n must be in the form 2^m-1 !
mod74

%% *** build laplacian matrix *********************************************
[ MLAP ] = laplacian74(n,nn);
disp('matrix ok')

%% *** initial function ***************************************************
L=6;
M=5;
[ ht_fI ,lape_fI   ]  = sol_exacte(x_fI,   y_fI,   z_fI   , L, M);
[ ht_fII ,lape_fII ]  = sol_exacte(x_fII,  y_fII,  z_fII  , L, M);
[ ht_fIII,lape_fIII]  = sol_exacte(x_fIII, y_fIII, z_fIII , L, M);
[ ht_fIV,lape_fIV  ]  = sol_exacte(x_fIV,  y_fIV,  z_fIV  , L, M);
[ ht_fV,lape_fV    ]  = sol_exacte(x_fV,   y_fV,   z_fV   , L, M);
[ ht_fVI,lape_fVI  ]  = sol_exacte(x_fVI,  y_fVI,  z_fVI  , L, M);

%% *** test the laplacian *************************************************
t1=cputime;
[ ht ] = resh( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn );
lap = MLAP*ht;
[ lap_fI, lap_fII, lap_fIII, lap_fIV, lap_fV, lap_fVI ] = deresh( lap, n, nn );
disp('step 1 : ok')
t_new=cputime-t1;
disp(['time (new method) : ' num2str(t_new)]);

t1=cputime;
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=gr74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn);
[lap2_fI,lap2_fII,lap2_fIII,lap2_fIV,lap2_fV,lap2_fVI]=div74(grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,n,nn);
disp('step 2 : ok')
t_old=cputime-t1;
disp(['time (old method) : ' num2str(t_old)]);

%% *** error **************************************************************
err_fI=lap_fI-lape_fI;
err_fII=lap_fII-lape_fII;
err_fIII=lap_fIII-lape_fIII;
err_fIV=lap_fIV-lape_fIV;
err_fV=lap_fV-lape_fV;
err_fVI=lap_fVI-lape_fVI;

%% *** curve **************************************************************
figure(1)
plot_cs7(n,nn,abs(err_fI), abs(err_fII), abs(err_fIII)...
    , abs(err_fIV), abs(err_fV), abs(err_fVI));
title('erreur with matricial laplacian')

figure(2)
plot_cs7(n,nn,abs(lap2_fI-lape_fI), abs(lap2_fII-lape_fII), abs(lap2_fIII-lape_fIII)...
    , abs(lap2_fIV-lape_fIV), abs(lap2_fV-lape_fV), abs(lap2_fVI-lape_fVI));
title('erreur with operator laplacian')

figure(3)
spy(MLAP)
title('matrix structure')

figure(4)
plot_cs7(n,nn,abs(lap_fI), abs(lap_fII), abs(lap_fIII)...
    , abs(lap_fIV), abs(lap_fV), abs(lap_fVI));
title('matricial laplacian')

figure(5)
plot_cs7(n,nn,abs(lap2_fI), abs(lap2_fII), abs(lap2_fIII),...
    abs(lap2_fIV), abs(lap2_fV), abs(lap2_fVI));
title('old laplacian')

figure(6)
plot_cs7(n,nn,abs(lape_fI), abs(lape_fII), abs(lape_fIII)...
    , abs(lape_fIV), abs(lape_fV), abs(lape_fVI));
title('exact laplacian')

figure(7)
plot_cs7(n,nn,abs(ht_fI), abs(ht_fII), abs(ht_fIII)...
    , abs(ht_fIV), abs(ht_fV), abs(ht_fVI));
title('norm of harmonic')

figure(8)
surf(y_fI./x_fI, z_fI./x_fI, abs(lap_fI-lape_fI))

fig_placier

%% *** evaluation de la constante *****************************************
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,lap_int]=...
    nrm74(lap_fI, lap_fII, lap_fIII, lap_fIV, lap_fV, lap_fVI,n,nn,'1');

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,ht_int]=...
    nrm74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn,'1');

eval=lap_int./ht_int;
eval_theo=L.*(L+1);
err=abs(eval-eval_theo)./abs(eval_theo)