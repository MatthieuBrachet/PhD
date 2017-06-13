clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global dxi
global opt_ftr scheme filtre

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';


test=2;
NN=2.^[2:7]-1;

for ii=1:length(NN)
    clc; n=NN(ii)
    mod101
    disp('mod74 : ok')
    
    hh(ii)=dxi;
    [fun_I    ,int] = fun_quad(x_fI  ,y_fI  ,z_fI   ,test);
    [fun_II   ,int] = fun_quad(x_fII ,y_fII ,z_fII  ,test);
    [fun_III  ,int] = fun_quad(x_fIII,y_fIII,z_fIII ,test);
    [fun_IV   ,int] = fun_quad(x_fIV ,y_fIV ,z_fIV  ,test);
    [fun_V    ,int] = fun_quad(x_fV  ,y_fV  ,z_fV   ,test);
    [fun_VI   ,int] = fun_quad(x_fVI ,y_fVI ,z_fVI  ,test);

    str='int';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    eint(ii)=abs((nrmI-int)/int);
    str='test';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    etest(ii)=abs((nrmI-int)/int);
    str='simpson';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg3]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    esimpson(ii)=abs((nrmI-int)/int);
end

[a11,b1]=polyfit(log(hh),log(eint),1);
[a12,b1]=polyfit(log(hh),log(etest),1);
[a13,b1]=polyfit(log(hh),log(esimpson),1);

clc
figure(1)
loglog(hh,eint,'-',hh,etest,'-',hh,esimpson,'-','LineWidth',2.0);
legend(['JPC rules              : ' num2str(a11(1))],['Trapezoidal rules : ' num2str(a12(1))],['Simpson rules      : ' num2str(a13(1))],'Location','SouthEast')
grid on

figure(2)
plot_cs11(n,nn,real(fun_I),real(fun_II),real(fun_III),real(fun_IV),real(fun_V),real(fun_VI))
%colormap colorcube
colorbar
title('real part')

figure(3)
plot_cs104(n,nn,imag(fun_I),imag(fun_II),imag(fun_III),imag(fun_IV),imag(fun_V),imag(fun_VI))
%colormap colorcube
colorbar
title('imag. part')

figure(4)
plot_cs104(n,nn,real(fun_I),real(fun_II),real(fun_III),real(fun_IV),real(fun_V),real(fun_VI))
%colormap colorcube
colorbar

fig_placier