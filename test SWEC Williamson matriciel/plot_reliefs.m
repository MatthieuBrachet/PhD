function []=plot_reliefs(n,nn)
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
global lambdac_mount tetac_mount

[hs_fI]   = relief(x_fI,y_fI,z_fI);
[hs_fII]  = relief(x_fII,y_fII,z_fII);
[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
[hs_fIV]  = relief(x_fIV,y_fIV,z_fIV);
[hs_fV]   = relief(x_fV,y_fV,z_fV);
[hs_fVI]  = relief(x_fVI,y_fVI,z_fVI);

hold on
%plot_cs100(n,nn,hs_fI,hs_fII,hs_fIII,hs_fIV,hs_fV,hs_fVI);

VThetaDeg = 0:.001:2*pi;
VTheta = VThetaDeg;
a=lambdac_mount;
b=tetac_mount;
R=pi/9;
XCercle = a + R * cos(VTheta);
YCercle = b + R * sin(VTheta);
ZCercle=5000*ones(size(XCercle));
plot3(XCercle, YCercle,ZCercle,'k-','Linewidth',2) 
view([0 0 10000])
hold off

