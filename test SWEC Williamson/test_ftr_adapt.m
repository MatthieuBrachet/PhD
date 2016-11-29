clc; clear all; close all

opt_ftr='redonnet10';
n=100;
h=1./n;
x=[h:h:1]';

opt_det='redonnet10';
opt_ftr='redonnet4';
[ detec ] = filtre74( n , opt_det );
[ ftr ] = filtre74( n , opt_ftr );
[ LAP_adap, B_adap, ftr1 ] = adaptative74( n, opt_ftr );

epsilon=10^-16;
rth=100;


E=[];
N=[];
for i=1:floor(n/2)
    y=cos(2*pi*i*x);
    
    uhf=y-detec*y;
    r=0.5*((B_adap*uhf).^2+((B_adap')*uhf).^2)./h.^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    uf=ftr*y;
    yf10=(1-sigma).*y+sigma.*uf;
    
    str=2;
    ee=norm(yf10,str)/norm(y,str);
    E=[E ee];
    
    
    N=[N i];
end

figure(1)
plot((2*pi/n)*N,E,'x-');
hold on


%opt_ftr='redonnet10';
n=100;
h=1./n;
x=[h:h:1]';

[ ftr ] = filtre74( n , opt_ftr );
[ ftr10 ] = filtre74( n , 'redonnet10' );

E=[];
E10=[];
N=[];
for i=1:floor(n/2)-1
    y=cos(2*pi*i*x);

    yf=ftr*y;
    
    str=2;
    ee=norm(yf,str)/norm(y,str);
    E=[E ee];
    
    yf10=ftr10*y;
    
    str=2;
    ee=norm(yf10,str)/norm(y,str);
    E10=[E10 ee];
    
    
    N=[N i];
end

figure(1)
plot((2*pi/n)*N,E,(2*pi/n)*N,E10)
legend('filtre adaptatif',opt_ftr,'redonnet 10')
