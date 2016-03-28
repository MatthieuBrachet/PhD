clc; clear all; close all;

A=[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16;15 17 18 19];
fun=A.^2+4;
betcr=A+0.5;;
[n1,n2]=size(A);

%% déconstruction en zigzag
bet=[];
betf=[];
betcr2=[];
for i=1:2:2*floor(n1/2)
   bet=[bet A(i,:) A(i+1,n2:-1:1)];
   betf=[betf fun(i,:) fun(i+1,n2:-1:1)];
   betcr2=[betcr2 betcr(i,:) betcr(i+1,n2:-1:1)];
end

ppspline=spline(bet,betf);
interp=ppval(ppspline,betcr2);

%% reconstruction
AA=[];
for i=1:2:n2
    AA=[AA; bet((i-1)*n1+1:i*n1); bet((i+1)*n1:-1:i*n1+1)];
end



% plot(bet,betf,'b-',betcr2,interp,'r-')
% legend('init','new')