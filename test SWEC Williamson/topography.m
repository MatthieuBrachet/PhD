function [ alt ] = topography( x,y,z )
load topo
rel=[topo topo];
[lambda,teta,rr]=cart2sph(x,y,z);
xx=pi/180*[1:1:720]-2*pi;
yy=pi/180*[1:180]-pi/2;
[X,Y]=meshgrid(xx,yy);
altp=interp2(X,Y,rel,lambda,teta);

sud=mean(topo(1,:));
nord=mean(topo(end,:));

[n1,n2]=size(x);
for i=1:n1
    for j=1:n2
        if abs(z(i,j)-rr)<10^-3
            alt(i,j)=nord;
        elseif abs(z(i,j)+rr)<10^-3
            alt(i,j)=sud;
        else
            alt(i,j)=altp(i,j);
        end
    end
end