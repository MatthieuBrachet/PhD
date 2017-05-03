function [ M ] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI )
%mass
global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI


[hs_fI] = relief(x_fI,y_fI,z_fI);
hstar_fI=ht_fI;%-hs_fI;

[hs_fII] = relief(x_fII,y_fII,z_fII);
hstar_fII=ht_fII;%-hs_fII;

[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
hstar_fIII=ht_fIII;%-hs_fIII;

[hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
hstar_fIV=ht_fIV;%-hs_fIV;

[hs_fV] = relief(x_fV,y_fV,z_fV);
hstar_fV=ht_fV;%-hs_fV;

[hs_fVI] = relief(x_fVI,y_fVI,z_fVI);
hstar_fVI=ht_fVI;%-hs_fVI;

nrm='cons_int';
%nrm='int';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,M]=nrm74(hstar_fI,hstar_fII,hstar_fIII,hstar_fIV,hstar_fV,hstar_fVI,n,nn,nrm);

end

