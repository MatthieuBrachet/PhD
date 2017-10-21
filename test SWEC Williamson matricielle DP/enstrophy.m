function [PE] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
% Potential enstrophy for SWE equation.
global nrm
global n nn
global omega alpha gp
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

[~, teta_fI, lambda_fI]=cart2sph(x_fI, y_fI, z_fI);
[~, teta_fII, lambda_fII]=cart2sph(x_fII, y_fII, z_fII);
[~, teta_fIII, lambda_fIII]=cart2sph(x_fIII, y_fIII, z_fIII);
[~, teta_fIV, lambda_fIV]=cart2sph(x_fIV, y_fIV, z_fIV);
[~, teta_fV, lambda_fV]=cart2sph(x_fV, y_fV, z_fV);
[~, teta_fVI, lambda_fVI]=cart2sph(x_fVI, y_fVI, z_fVI);

cor_fI=2.*omega.*(-cos(lambda_fI).*cos(teta_fI).*sin(alpha)+sin(teta_fI).*cos(alpha));
cor_fII=2.*omega.*(-cos(lambda_fII).*cos(teta_fII).*sin(alpha)+sin(teta_fII).*cos(alpha));
cor_fIII=2.*omega.*(-cos(lambda_fIII).*cos(teta_fIII).*sin(alpha)+sin(teta_fIII).*cos(alpha));
cor_fIV=2.*omega.*(-cos(lambda_fIV).*cos(teta_fIV).*sin(alpha)+sin(teta_fIV).*cos(alpha));
cor_fV=2.*omega.*(-cos(lambda_fV).*cos(teta_fV).*sin(alpha)+sin(teta_fV).*cos(alpha));
cor_fVI=2.*omega.*(-cos(lambda_fVI).*cos(teta_fVI).*sin(alpha)+sin(teta_fVI).*cos(alpha));

for i=1:n1
    for j=1:n2
        xi_fI(i,j)=0.5*((vort_fI(i,j)+cor_fI(i,j)).^2)./(gp*hstar_fI(i,j));
        xi_fII(i,j)=0.5*((vort_fII(i,j)+cor_fII(i,j)).^2)./(gp*hstar_fII(i,j));
        xi_fIII(i,j)=0.5*((vort_fIII(i,j)+cor_fIII(i,j)).^2)./(gp*hstar_fIII(i,j));
        xi_fIV(i,j)=0.5*((vort_fIV(i,j)+cor_fIV(i,j)).^2)./(gp*hstar_fIV(i,j));
        xi_fV(i,j)=0.5*((vort_fV(i,j)+cor_fV(i,j)).^2)./(gp*hstar_fV(i,j));
        xi_fVI(i,j)=0.5*((vort_fVI(i,j)+cor_fVI(i,j)).^2)./(gp*hstar_fVI(i,j));
    end
end

[PE_fI,PE_fII,PE_fIII,PE_fIV,PE_fV,PE_fVI,PE]=...
    nrm101(xi_fI,xi_fII,xi_fIII,xi_fIV,xi_fV,xi_fVI,n,nn,nrm);


end

