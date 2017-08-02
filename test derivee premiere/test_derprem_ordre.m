clc; clear all; close all
NN=10:10:100;
err=[];hh=[];
for i=1:length(NN)
    N=NN(i);
    x=linspace(0,1,N)';
    h=x(2)-x(1);

    y=cos(x)+x.^2-3;
    dye=-sin(x)+2*x;

    P=diag(4/6*ones(N,1))+diag(1/6*ones(N-1,1),1)+diag(1/6*ones(N-1,1),-1);
    P=sparse(P);

    a=-103/72; b=91/36; c=-7/4; d=29/36; e=-11/72;
    Q=(diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/2;
    Q(1,1)=a;
    Q(1,2)=b;
    Q(1,3)=c;
    Q(1,4)=d;
    Q(1,5)=e;
    Q(end,end)=-a;
    Q(end,end-1)=-b;
    Q(end,end-2)=-c;
    Q(end,end-3)=-d;
    Q(end,end-4)=-e;
    Q=sparse(Q)/h;

    w=Q*y;
    dy=P\w;

    err=[err norm(dy-dye)];
    hh=[hh h];
end

figure(1)
loglog(hh,err,hh,hh.^4)