function [Ec] = kineticenergy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI

[hs_fI] = relief(x_fI,y_fI,z_fI);
hstar_fI=ht_fI-hs_fI;

[hs_fII] = relief(x_fII,y_fII,z_fII);
hstar_fII=ht_fII-hs_fII;

[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
hstar_fIII=ht_fIII-hs_fIII;

[hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
hstar_fIV=ht_fIV-hs_fIV;

[hs_fV] = relief(x_fV,y_fV,z_fV);
hstar_fV=ht_fV-hs_fV;

[hs_fVI] = relief(x_fVI,y_fVI,z_fVI);
hstar_fVI=ht_fVI-hs_fVI;

[n1,n2]=size(x_fI);
for i=1:n1
    for j=1:n2
        norm_I(i,j)=dot(hstar_fI(i,j)*vt_fI(i,j,:),vt_fI(i,j,:));
        norm_II(i,j)=dot(hstar_fII(i,j)*vt_fII(i,j,:),vt_fII(i,j,:));
        norm_III(i,j)=dot(hstar_fIII(i,j)*vt_fIII(i,j,:),vt_fIII(i,j,:));
        norm_IV(i,j)=dot(hstar_fIV(i,j)*vt_fIV(i,j,:),vt_fIV(i,j,:));
        norm_V(i,j)=dot(hstar_fV(i,j)*vt_fV(i,j,:),vt_fV(i,j,:));
        norm_VI(i,j)=dot(hstar_fVI(i,j)*vt_fVI(i,j,:),vt_fVI(i,j,:));
    end
end
nrm='cor_int';
[~,~,~,~,~,~,totalu]=...
    nrm101(norm_I,norm_II,norm_III,norm_IV,norm_V,norm_VI,n,nn,nrm);

Ec=0.5*totalu;


end