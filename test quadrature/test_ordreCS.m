clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global dxi corner
global opt_ftr scheme filtre

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

corner = -100;

test=7;
NN=2.^[3:6]-1;

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

    str='uniforme';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg0]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    euni(ii)=abs((nrmg0-int));%/int);
    str='int';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    eint(ii)=abs((nrmg1-int));%/int);
    str='trapezes';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    erect(ii)=abs((nrmg2-int));%/int);
    str='simpson';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg3]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    esimpson(ii)=abs((nrmg3-int));%/int);
    str='test';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg4]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
    etest(ii)=abs((nrmg4-int));%/int);
end

[a10,b1]=polyfit(log(hh),log(euni),1);
[a11,b1]=polyfit(log(hh),log(eint),1);
[a12,b1]=polyfit(log(hh),log(erect),1);
[a13,b1]=polyfit(log(hh),log(esimpson),1);
[a14,b1]=polyfit(log(hh),log(etest),1);

clc;
figure(1)
loglog(hh,euni,'-',hh,erect,'-',hh,etest,'-',hh,eint,'-',hh,esimpson,'-','LineWidth',2.0);
grid on
xlabel('\Delta \xi')
ylabel('Error on quadrature')
legend(['Uniforme rules           : ' num2str(a10(1))],['Trapezoidal rules       : ' num2str(a12(1))],['Q_{\alpha} rules with \alpha=' num2str(corner) '     : ' num2str(a14(1))],['Q_{\alpha} rules with \alpha=1/3  : ' num2str(a11(1))],['Simpson rules           : ' num2str(a13(1))],'Location','SouthEast')

figure(2)
plot_cs11(n,nn,real(fun_I),real(fun_II),real(fun_III),real(fun_IV),real(fun_V),real(fun_VI))
colorbar
title('real part')

figure(3)
plot_cs104(n,nn,imag(fun_I),imag(fun_II),imag(fun_III),imag(fun_IV),imag(fun_V),imag(fun_VI))
colorbar
title('imag. part')

figure(4)
plot_cs104(n,nn,real(fun_I),real(fun_II),real(fun_III),real(fun_IV),real(fun_V),real(fun_VI))
colorbar
title('real part')

fig_placier
