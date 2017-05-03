function [qt_x,qt_y,qt_z] = qt_mvt(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
global n nn
%% x-comp
qt_x_fI=ht_fI.*vt_fI(:,:,1);
qt_x_fII=ht_fII.*vt_fII(:,:,1);
qt_x_fIII=ht_fIII.*vt_fIII(:,:,1);
qt_x_fIV=ht_fIV.*vt_fIV(:,:,1);
qt_x_fV=ht_fV.*vt_fV(:,:,1);
qt_x_fVI=ht_fVI.*vt_fVI(:,:,1);

nrm='cor_int';
[~,~,~,~,~,~,qt_x]=...
    nrm101(qt_x_fI,qt_x_fII,qt_x_fIII,qt_x_fIV,qt_x_fV,qt_x_fVI,n,nn,nrm);

%% y-comp
qt_y_fI=ht_fI.*vt_fI(:,:,2);
qt_y_fII=ht_fII.*vt_fII(:,:,2);
qt_y_fIII=ht_fIII.*vt_fIII(:,:,2);
qt_y_fIV=ht_fIV.*vt_fIV(:,:,2);
qt_y_fV=ht_fV.*vt_fV(:,:,2);
qt_y_fVI=ht_fVI.*vt_fVI(:,:,2);

nrm='cor_int';
[~,~,~,~,~,~,qt_y]=...
    nrm101(qt_y_fI,qt_y_fII,qt_y_fIII,qt_y_fIV,qt_y_fV,qt_y_fVI,n,nn,nrm);

%% z-comp
qt_z_fI=ht_fI.*vt_fI(:,:,3);
qt_z_fII=ht_fII.*vt_fII(:,:,3);
qt_z_fIII=ht_fIII.*vt_fIII(:,:,3);
qt_z_fIV=ht_fIV.*vt_fIV(:,:,3);
qt_z_fV=ht_fV.*vt_fV(:,:,3);
qt_z_fVI=ht_fVI.*vt_fVI(:,:,3);

nrm='cor_int';
[~,~,~,~,~,~,qt_z]=...
    nrm101(qt_z_fI,qt_z_fII,qt_z_fIII,qt_z_fIV,qt_z_fV,qt_z_fVI,n,nn,nrm);

end