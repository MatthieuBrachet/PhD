clc; clear all; close all;

a=0;
b=1;

epsilon=10^-4;
itermax=100;

n=0;
e=1;
E=[];
ERR=[];
ITER=[];

I=(fun(a)+fun(b))/2;
while e>epsilon & n<itermax
    n=n+1;
    
    S=0;
    for ppp=0:n-1
        S=S+fun(a+(2*ppp+1)*(b-a)/(2^n));
    end
    In=I+S;
    
    ITER=[ITER n];
    
    e=abs(In-I)/abs(I);
    E=[E e];
    
    I=In;
    INTEG(n)=I*(b-a)/(2^n);
    ERREUR(n)=abs(INTEG(n)-1/3);
end

figure(1)
loglog(ITER,E,ITER,ERREUR)
legend('relativ error');
    

figure(2)
plot(ITER,INTEG)