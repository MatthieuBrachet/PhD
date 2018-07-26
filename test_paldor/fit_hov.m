function [ Cn ] = fit_hov( X,T,U )
[p,q]=size(X);
for j=1:1
    xx=T(:,q);
    pts=U(:,q);
    f = fit(xx,pts,'fourier8');
    cn(j)=f.w;
end
Cn=mean(cn);