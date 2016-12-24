clc; clear all; close all

% f=cox(x)
x=0:0.01:pi;
y=sin(x);
plot(x,y)

nnn=10000;
a=0;
b=pi;
fa=sin(a);
fb=sin(b);
h=(b-a)/nnn;

s1=0;
for kk=1:nnn/2-1
    x=a+(2*kk)*h;
    fx=sin(x);
    s1=s1+fx;
end

s2=0;
for kk=1:nnn/2
    x=a+(2*kk-1)*h;
    fx=sin(x);
    s2=s2+fx;
end

I=(h/3)*(fa+2*s1+4*s2+fb)