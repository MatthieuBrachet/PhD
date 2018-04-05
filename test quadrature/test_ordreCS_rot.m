clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global dxi corner
global opt_ftr scheme filtre
global ang1 ang2

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
corner = 1;

test=0;
NN=2.^[3:8]-1;

for ii=1:length(NN)
    n=NN(ii);
    mod101
    
    eint(ii)=0;
    erect(ii)=0;
    esimpson(ii)=0;
    etest(ii)=0;
    for iter = 1:1000
        clc; [n iter]
        ang1=2*pi*rand(1);
        ang2=2*pi*rand(1);
    
        hh(ii)=dxi;
        [fun_I    ,int] = fun_quad2(x_fI  ,y_fI  ,z_fI   ,test);
        [fun_II   ,int] = fun_quad2(x_fII ,y_fII ,z_fII  ,test);
        [fun_III  ,int] = fun_quad2(x_fIII,y_fIII,z_fIII ,test);
        [fun_IV   ,int] = fun_quad2(x_fIV ,y_fIV ,z_fIV  ,test);
        [fun_V    ,int] = fun_quad2(x_fV  ,y_fV  ,z_fV   ,test);
        [fun_VI   ,int] = fun_quad2(x_fVI ,y_fVI ,z_fVI  ,test);

        str='int';
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
        eint(ii)=max(eint(ii),abs((nrmg1-int))/int);
        str='trapezes';
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
        erect(ii)=max(erect(ii),abs((nrmg2-int))/int);
        str='simpson';
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg3]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
        esimpson(ii)=max(esimpson(ii),abs((nrmg3-int))/int);
        str='test';
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg4]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
        etest(ii)=max(etest(ii),abs((nrmg4-int))/int);
    end
end


[a11,b11]=polyfit(log10(hh),log10(eint),1);
[a12,b12]=polyfit(log10(hh),log10(erect),1);
[a13,b13]=polyfit(log10(hh),log10(esimpson),1);
[a14,b14]=polyfit(log10(hh),log10(etest),1);




figure(1)
hold on
hl1=plot(log10(hh),log10(eint),'kx','LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'LineWidth',2.0);
hlf1=plot(log10(hh),a11(1)*log10(hh)+a11(2),'b-','LineWidth',2.0);

hl2=plot(log10(hh),log10(erect),'kx','LineWidth',2.0);
set(hl2,'MarkerSize',10);
set(hl2,'LineWidth',2.0);
hlf2=plot(log10(hh),a12(1)*log10(hh)+a12(2),'r-.','LineWidth',2.0);

hl3=plot(log10(hh),log10(esimpson),'kx','LineWidth',2.0);
set(hl3,'MarkerSize',10);
set(hl2,'LineWidth',2.0);
hlf3=plot(log10(hh),a13(1)*log10(hh)+a13(2),'m--','LineWidth',2.0);

hl4=plot(log10(hh),log10(etest),'kx','LineWidth',2.0);
set(hl4,'MarkerSize',10);
set(hl4,'LineWidth',2.0);
hlf4=plot(log10(hh),a14(1)*log10(hh)+a14(2),'g:','LineWidth',2.0);
hold off

txt1=['Q_{\alpha} avec \alpha=1/3  - pente : ' num2str(a11(1))];
txt2=['Trapèzes           - pente  : ' num2str(a12(1))];
txt3=['Simpson            - pente  : ' num2str(a13(1))];
txt4=['Q_{\alpha} avec \alpha=' num2str(corner) '     - pente  : ' num2str(a14(1))];


legend([hlf4,hlf2,hlf1,hlf3],{txt4,txt2,txt1,txt3},'Location','SouthEast')
grid on
xlabel('Log_{10}(\Delta)')
ylabel('Log_{10}(Erreur)')
axis([-2.5 -.5 -11 -1])