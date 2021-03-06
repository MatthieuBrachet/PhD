function [y] = fun10(x, teta0, teta1)
% function extract from J. Galewsky, R.K. Scott and L.M. Polvani article
% (an initial-value problem for testing numerical models of the global
% shallow water equations).
[n1,n2]=size(x);
y=zeros(n1,n2);
en=exp(-4/(teta0-teta1)^2);
for i=1:n1
    for j=1:n2
        if x(i,j)<=teta0
            y(i,j)=0;
        elseif x(i,j)>=teta1
            y(i,j)=0;
        else
            y(i,j)=exp(1/((x(i,j)-teta0)*(x(i,j)-teta1)))./en;
        end
    end
end