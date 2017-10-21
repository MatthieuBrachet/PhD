function [Epp] = Potentialenergy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
global gp

[hs_fI] = relief(x_fI,y_fI,z_fI);
[hs_fII] = relief(x_fII,y_fII,z_fII);
[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
[hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
[hs_fV] = relief(x_fV,y_fV,z_fV);
[hs_fVI] = relief(x_fVI,y_fVI,z_fVI);

nrm='cor_int';
[~,~,~,~,~,~,totalh]=...
    nrm101(ht_fI.^2-hs_fI.^2,ht_fII.^2-hs_fII.^2,ht_fIII.^2-hs_fIII.^2,ht_fIV.^2-hs_fIV.^2,ht_fV.^2-hs_fV.^2,ht_fVI.^2-hs_fVI.^2,n,nn,nrm);

Epp=0.5*gp*totalh;


end