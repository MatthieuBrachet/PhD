clc; clear all; close all

opt_ftr='redonnet10';
n=100;
h=1./n;
x=[h:h:1]';

[ ftr ] = filtre74( n , opt_ftr );
[ ftr10 ] = filtre74( n , 'redonnet6' );

E=[];
E10=[];
N=[];
for i=1:floor(n/2)-1
    y=cos(2*pi*i*x);

    yf=ftr*y;
    
    str=2;
    ee=norm(y-yf,str)/norm(y,str);
    E=[E ee];
    
    yf10=ftr10*y;
    
    str=2;
    ee=norm(y-yf10,str)/norm(y,str);
    E10=[E10 ee];
    
    
    N=[N i];
end

figure(2)
semilogy((2*pi/n)*N,E,(2*pi/n)*N,E10)
legend('bogey 6','redonnet 6')
