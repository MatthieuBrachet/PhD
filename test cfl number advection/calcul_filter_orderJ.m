clc; clear all; close all;
format long
% calcul d'un filtre symétrique d'ordre 2*J avec J>0.
J=21;

A=zeros(J+1,J+1);
NN=[0:J];
A(1,:)=1;
A(2,:)=(-1).^NN;
for i=3:J+1
    for j=2:J+1
        A(i,j)=(j-1).^(2*(i-2));
    end
end
b=zeros(J+1,1);
b(1)=1;
%[a,FLAG,RELRES,ITER]=gmres(A,b);
P=tril(A);
Ap=P\A;
bp=P\b;
a=Ap\bp;
err=norm(A*a-b)
sum(a)

figure(2)
plot(a)

teta=linspace(0,pi,10000);
ampli=zeros(size(teta));
for k=1:length(a)
    ampli=ampli+a(k).*cos((k-1)*teta);
end

figure(1)
plot(teta,ampli)

max(abs(ampli))
