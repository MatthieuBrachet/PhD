function [h,u,v] = fun100(X,Y)
h0=1;
r=sqrt((X-.5).^2+(Y-.5).^2);
h=h0.*exp(-(r.^2./(.01)));
u=zeros(size(h));
v=zeros(size(h));
h=reshape(h,[],1);
u=reshape(u,[],1);
v=reshape(v,[],1);
end