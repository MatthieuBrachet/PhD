clc; clear all; close all

opt_ftr='bogey6';
NN=10:10:200;

E=[];
E6=[];
H=[];
for i=1:length(NN)
    clc; nn=NN(i)
    h=1./nn;
    H=[H h];
    
    x=[h:h:1]';
    y=cos(2*pi*x);
    
    [ ftr ] = filtre74( nn , opt_ftr );
    yf=ftr*y;
    ee=norm(y-yf,2);
    E=[E ee];
    
    
    [ ftr6 ] = filtre74( nn , 'redonnet6' );
    yf6=ftr6*y;
    ee=norm(y-yf6,2);
    E6=[E6 ee];
    
end

figure(1)
loglog(H,E,H,H.^6,H,E6)
legend('numerical','theorical','referential order 6')
grid on

P=polyfit(log(H),log(E),1);
disp(P(1))

figure(2)
plot(x,y,'o-',x,yf,'x-')
legend('pert. func.','filtered func')
