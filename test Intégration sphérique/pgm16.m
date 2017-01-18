%%% Extrapolation de Richardson sur la formule sans les epsilon

clear all;
global n nn;
global radius;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

int1=[];
int2=[];
int3=[];
Nmax=120;
tmp=size(4:2:Nmax);
s=tmp(2);
% 216*pi/35
% 6.6961822200736179523

for N=4:2:Nmax,

    N
    
    make_cs_grid(N);
    funfI=fun_f2(x_fI,y_fI,z_fI);
    funfII=fun_f2(x_fII,y_fII,z_fII);
    funfIII=fun_f2(x_fIII,y_fIII,z_fIII);
    funfIV=fun_f2(x_fIV,y_fIV,z_fIV);
    funfV=fun_f2(x_fV,y_fV,z_fV);
    funfVI=fun_f2(x_fVI,y_fVI,z_fVI);
    weights=dxi*deta*dga;
    [tmp1,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
   
    int3=[int3;tmp1];
    
    N=2*N;
    make_cs_grid(N);
    funfI=fun_f2(x_fI,y_fI,z_fI);
    funfII=fun_f2(x_fII,y_fII,z_fII);
    funfIII=fun_f2(x_fIII,y_fIII,z_fIII);
    funfIV=fun_f2(x_fIV,y_fIV,z_fIV);
    funfV=fun_f2(x_fV,y_fV,z_fV);
    funfVI=fun_f2(x_fVI,y_fVI,z_fVI);
    weights=dxi*deta*dga;
    [tmp2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
   
    int1=[int1;abs((2^4*tmp2-tmp1)/(2^4-1))];
    
end

err1=abs(int1-6.6961822200736179523);
err3=abs(int3-6.6961822200736179523);
polyfit(log(20:2:Nmax)',log(err1(9:end)),1)


    
for i=1:(s-1)/2,
    int2=[int2;abs((2^6*int1(2*i+1)-int1(i))/(2^6-1))];
end

err2=abs(int2-6.6961822200736179523);
polyfit(log(20:2:Nmax/2-5*2)',log(err2(9:end-5)),1)

    
p1=plot(log(4:2:Nmax),log(err1),'-s');
grid minor; grid on;
hold on;
p2=plot(log(4:2:Nmax/2),log(err2),'-d');
hold on;
x=[4;8;16;32;64]; 
y=[10^(-3);10^(-6);10^(-10);10^(-14);10^(-18)];
%y=[10^(-4);10^(-15);10^(-18);10^(-18);10^(-18)];
p3=plot(log(x),log(y),'-');
p3.LineWidth=1.5;
hold on;
p4=plot(log(4:2:Nmax),log(err3),'-');
legend('Richardson 1','Richardson 2','Epsilons','De base','Location','southwest');
title('Log erreurs intégration pour f_2');
xlabel('log(N)'); ylabel('log(erreur)');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    