function [ I ] = simpson(f,a,b,n)
h=(b-a)/n;
s1=0;
s2=0;
for i=1:floor(n/2)-1
    aa=a+2*i*h;
    s1=s1+f(aa);
end
for i=1:floor(n/2)
    bb=a+(2*i-1)*h;
    s2=s2+f(bb);
end
I=h/3*(f(a)+2*s1+4*s2+f(b));
end