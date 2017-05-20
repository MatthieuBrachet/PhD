clc; clear all; close all;

NN=[1:2:10 20:100:1000];
E=[]; H=[]; INT=[];
for i=1:length(NN);
    n=NN(i);
    h=1/(n+1);
    x=0:h:1;
    [X,Y]=meshgrid(x,x);
    F=(sin(X)+X).*Y;

    wei=ones(size(X));
    border=1/2;
    corner=1/80;
    
    wei(1,1:end)=border;
    wei(end,1:end)=border;
    wei(1:end,1)=border;
    wei(1:end,end)=border;
    wei(1,1)=corner;
    wei(1,end)=corner;
    wei(end,1)=corner;
    wei(end,end)=corner;

    int=h^2.*sum(sum(wei.*F));
    inte=.5*(1.5-cos(1));
    E=[E int-inte];
    H=[H h];
    INT=[INT int];
end

figure(1)
surf(X,Y,F)
shading interp

figure(2)
loglog(H,E)

figure(3)
loglog(INT)

p=polyfit(log(H),log(abs(E)),1);
disp(['ordre estimé : ' num2str(p(1))])

