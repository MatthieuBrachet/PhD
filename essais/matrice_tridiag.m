clc; clear all; close all;


alfa=0.3;
beta=4;
T1=[];
T2=[];
NUM=[];

NN=10:50:500;
for pp=1:10
    pp
    NUM=[NUM NN];
    for i=1:length(NN)
        n=NN(i);
    
        A=eye(n,n)+alfa*(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1));A=sparse(A);
        B=eye(n,n)+beta*(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1));B=sparse(B);
        C=A*B;
    
        sm=ones(n,1);
    
        tic;
        w=A\sm;
        x=B\w;
        time1=toc;
        T1=[T1 time1];
    
    
        tic;
        x=C\sm;
        time2=toc;
        T2=[T2 time2];
    end
end

figure(1)
plot(NUM,T1,'bx',NUM,T2,'rx')
legend('systemes separes','produit')