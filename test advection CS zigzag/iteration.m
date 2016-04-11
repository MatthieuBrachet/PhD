function [funfInew, funfIInew, funfIIInew, funfIVnew, funfVnew, funfVInew] = iteration(funfI, funfII,funfIII,funfIV,funfV,funfVI,ddt,time)

global n nn


[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_1b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

funfI=funftI;funfII=funftII;funfIII=funftIII;
funfIV=funftIV;funfV=funftV;funfVI=funftVI;

%%   CALCUL RK4
%    iterations

%% CALCUL KK0
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr_1b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time);

kk0_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kk0_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kk0_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kk0_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kk0_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kk0_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 


%% CALCUL KK1

fun1_fI=funfI+0.5*ddt*kk0_I;
fun1_fII=funfII+0.5*ddt*kk0_II;
fun1_fIII=funfIII+0.5*ddt*kk0_III;
fun1_fIV=funfIV+0.5*ddt*kk0_IV;
fun1_fV=funfV+0.5*ddt*kk0_V;
fun1_fVI=funfVI+0.5*ddt*kk0_VI;


[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr_1b(fun1_fI,fun1_fII,fun1_fIII,fun1_fIV,fun1_fV,fun1_fVI,n,nn);

[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt/2);


kk1_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kk1_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kk1_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kk1_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kk1_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kk1_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 


%% CALCUL KK2
fun2_fI=funfI+0.5*ddt*kk1_I;
fun2_fII=funfII+0.5*ddt*kk1_II;
fun2_fIII=funfIII+0.5*ddt*kk1_III;
fun2_fIV=funfIV+0.5*ddt*kk1_IV;
fun2_fV=funfV+0.5*ddt*kk1_V;
fun2_fVI=funfVI+0.5*ddt*kk1_VI;

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr_1b(fun2_fI,fun2_fII,fun2_fIII,fun2_fIV,fun2_fV,fun2_fVI,n,nn);

[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt/2);

kk2_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kk2_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kk2_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kk2_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kk2_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kk2_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 


%% CALCUL KK3
fun3_fI=funfI+ddt*kk2_I;
fun3_fII=funfII+ddt*kk2_II;
fun3_fIII=funfIII+ddt*kk2_III;
fun3_fIV=funfIV+ddt*kk2_IV;
fun3_fV=funfV+ddt*kk2_V;
fun3_fVI=funfVI+ddt*kk2_VI;


 [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
     gr_1b(fun3_fI,fun3_fII,fun3_fIII,fun3_fIV,fun3_fV,fun3_fVI,n,nn);

 [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt);
 
kk3_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kk3_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kk3_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kk3_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kk3_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kk3_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 


%% ASSEMBLAGE RK4
funfInew(1:nn,1:nn)=funfI(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_I(1:nn,1:nn)+(1/3)*kk1_I(1:nn,1:nn)+(1/3)*kk2_I(1:nn,1:nn)+(1/6)*kk3_I(1:nn,1:nn));
funfIInew(1:nn,1:nn)=funfII(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_II(1:nn,1:nn)+(1/3)*kk1_II(1:nn,1:nn)+(1/3)*kk2_II(1:nn,1:nn)+(1/6)*kk3_II(1:nn,1:nn));
funfIIInew(1:nn,1:nn)=funfIII(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_III(1:nn,1:nn)+(1/3)*kk1_III(1:nn,1:nn)+(1/3)*kk2_III(1:nn,1:nn)+(1/6)*kk3_III(1:nn,1:nn));
funfIVnew(1:nn,1:nn)=funfIV(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_IV(1:nn,1:nn)+(1/3)*kk1_IV(1:nn,1:nn)+(1/3)*kk2_IV(1:nn,1:nn)+(1/6)*kk3_IV(1:nn,1:nn));
funfVnew(1:nn,1:nn)=funfV(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_V(1:nn,1:nn)+(1/3)*kk1_V(1:nn,1:nn)+(1/3)*kk2_V(1:nn,1:nn)+(1/6)*kk3_V(1:nn,1:nn));
funfVInew(1:nn,1:nn)=funfVI(1:nn,1:nn)+ddt*...
    ((1/6)*kk0_VI(1:nn,1:nn)+(1/3)*kk1_VI(1:nn,1:nn)+(1/3)*kk2_VI(1:nn,1:nn)+(1/6)*kk3_VI(1:nn,1:nn));
%
[funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew]=...
    ds_1b(funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew,n,nn);


end

