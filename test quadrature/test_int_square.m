clc; clear all; close all;

NN=[2:2:10 20:100:2000]+1;
E=[]; H=[]; INT=[];
for i=1:length(NN);
    n=NN(i);
    h=1/(n+1);
    x=0:h:1;
    [X,Y]=meshgrid(x,x);
    
%     F=(cos(X)+X).*Y;
%     inte=.5*(.5+sin(1));
    
    pp=2;
    qq=3;
    F=X.^pp.*Y.^qq;
    inte=1/((pp+1)*(qq+1));

    wei=ones(size(X));
    border=1/2;
    corner=1/4;
    
    wei(1,1:end)=border;
    wei(end,1:end)=border;
    wei(1:end,1)=border;
    wei(1:end,end)=border;
    wei(1,1)=corner;
    wei(1,end)=corner;
    wei(end,1)=corner;
    wei(end,end)=corner;
    
    int=h^2.*sum(sum(wei.*F));

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
plot(INT-inte)

fig_placier

p=polyfit(log(H),log(abs(E)),1);
disp(['ordre estimé : ' num2str(p(1))])

