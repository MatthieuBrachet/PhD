function [ alt ] = topography( x,y,z )
load topo
rel=[topo topo];
[lambda,teta,rr]=cart2sph(x,y,z);
xx=pi/180*[1:1:720]-2*pi;
yy=pi/180*[1:180]-pi/2;
[X,Y]=meshgrid(xx,yy);
altp=interp2(X,Y,rel,lambda,teta);
sud=mean(topo(1,:));
altp(isnan(altp))=0;
alt=altp.*(1-(z==-rr))+sud.*(z==-rr);
alt=alt.*(alt>0);