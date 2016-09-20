clc; clear all; close all;
format long;

global n
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global alpha u0 radius

test=0;

alpha=pi/3;
video = 'no';
sauvegarde = 1;
opt_ftr=10;
scheme='compact4';
NN=10:20:100;

EQ1=[]; EQ2=[];
for i=1:length(NN)
    clc
    n=NN(i)
    mod74
    u0=2*pi*radius/(12*24*60*60);
    t=0;
    [ht_fI,vt_fI] = sol_exacte(x_fI,y_fI,z_fI,t);
    [ht_fII,vt_fII] = sol_exacte(x_fII,y_fII,z_fII,t);
    [ht_fIII,vt_fIII] = sol_exacte(x_fIII,y_fIII,z_fIII,t);
    [ht_fIV,vt_fIV] = sol_exacte(x_fIV,y_fIV,z_fIV,t);
    [ht_fV,vt_fV] = sol_exacte(x_fV,y_fV,z_fV,t);
    [ht_fVI,vt_fVI] = sol_exacte(x_fVI,y_fVI,z_fVI,t);


    [eq1_fI, eq1_fII, eq1_fIII, eq1_fIV, eq1_fV, eq1_fVI]=eq_moment74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
    [eq2_fI, eq2_fII, eq2_fIII, eq2_fIV, eq2_fV, eq2_fVI]=eq_cons74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);

    e_eq1=max(max(max(abs([eq1_fI, eq1_fII, eq1_fIII, eq1_fIV, eq1_fV, eq1_fVI]))));
    EQ1=[EQ1 e_eq1];
    e_eq2=max(max(abs([eq2_fI, eq2_fII, eq2_fIII, eq2_fIV, eq2_fV, eq2_fVI])));
    EQ2=[EQ2 e_eq2];
end

figure(1)
loglog(1./NN,EQ1,1./NN,EQ2,1./NN, 1./(NN.^4))
legend('eq moment','eq. cons.','h^4')
