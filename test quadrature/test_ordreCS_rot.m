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

test=2;
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

[a11,b1]=polyfit(log(hh),log(eint),1);
[a12,b1]=polyfit(log(hh),log(erect),1);
[a13,b1]=polyfit(log(hh),log(esimpson),1);
[a14,b1]=polyfit(log(hh),log(etest),1);

figure(1)
loglog(hh,erect,'-',hh,etest,'-',hh,eint,'-',hh,esimpson,'-','LineWidth',2.0);
legend(['Trapezoidal rules       : ' num2str(a12(1))],['Q_{\alpha} rules with \alpha=' num2str(corner) '     : ' num2str(a14(1))],['Q_{\alpha} rules with \alpha=1/3  : ' num2str(a11(1))],['Simpson rules            : ' num2str(a13(1))],'Location','SouthEast')
grid on
xlabel('\Delta \xi')
ylabel('Error on quadrature')
axis([hh(end)*.95 hh(1)*1.05 10^-11 10^0])
