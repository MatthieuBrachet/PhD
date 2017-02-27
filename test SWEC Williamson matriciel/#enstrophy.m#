function [PE] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
% Potential enstrophy
global n nn
global omega
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

[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
    vort101(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

[n1,n2]=size(x_fI);

[~, teta_fI, ~]=cart2sph(x_fI, y_fI, z_fI);
[~, teta_fII, ~]=cart2sph(x_fII, y_fII, z_fII);
[~, teta_fIII, ~]=cart2sph(x_fIII, y_fIII, z_fIII);
[~, teta_fIV, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
[~, teta_fV, ~]=cart2sph(x_fV, y_fV, z_fV);
[~, teta_fVI, ~]=cart2sph(x_fVI, y_fVI, z_fVI);

cor_fI=2.*omega.*sin(teta_fI);
cor_fII=2.*omega.*sin(teta_fII);
cor_fIII=2.*omega.*sin(teta_fIII);
cor_fIV=2.*omega.*sin(teta_fIV);
cor_fV=2.*omega.*sin(teta_fV);
cor_fVI=2.*omega.*sin(teta_fVI);

for i=1:n1
    for j=1:n2
        xi_fI(i,j)=0.5*((vort_fI(i,j)+cor_fI(i,j)).^2)./hstar_fI(i,j);
        xi_fII(i,j)=0.5*((vort_fII(i,j)+cor_fII(i,j)).^2)./hstar_fII(i,j);
        xi_fIII(i,j)=0.5*((vort_fIII(i,j)+cor_fIII(i,j)).^2)./hstar_fIII(i,j);
        xi_fIV(i,j)=0.5*((vort_fIV(i,j)+cor_fIV(i,j)).^2)./hstar_fIV(i,j);
        xi_fV(i,j)=0.5*((vort_fV(i,j)+cor_fV(i,j)).^2)./hstar_fV(i,j);
        xi_fVI(i,j)=0.5*((vort_fVI(i,j)+cor_fVI(i,j)).^2)./hstar_fVI(i,j);
    end
end

nrm='cor_int';
%nrm='int';
[~,~,~,~,~,~,PE]=...
    nrm101(xi_fI,xi_fII,xi_fIII,xi_fIV,xi_fV,xi_fVI,n,nn,nrm);


end

