function [ E ] = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn)
global gp hp

norm_I=vt_fI(:,:,1).^2+vt_fI(:,:,2).^2+vt_fI(:,:,3).^2;
norm_II=vt_fII(:,:,1).^2+vt_fII(:,:,2).^2+vt_fII(:,:,3).^2;
norm_III=vt_fIII(:,:,1).^2+vt_fIII(:,:,2).^2+vt_fIII(:,:,3).^2;
norm_IV=vt_fIV(:,:,1).^2+vt_fIV(:,:,2).^2+vt_fIV(:,:,3).^2;
norm_V=vt_fV(:,:,1).^2+vt_fV(:,:,2).^2+vt_fV(:,:,3).^2;
norm_VI=vt_fVI(:,:,1).^2+vt_fVI(:,:,2).^2+vt_fVI(:,:,3).^2;


[aaa,aaa,aaa,aaa,aaa,aaa,totalu]=...
    nrm101(norm_I,norm_II,norm_III,norm_IV,norm_V,norm_VI,n,nn,'cor_int');

[aaa,aaa,aaa,aaa,aaa,aaa,totalh]=...
    nrm101(ht_fI.^2,ht_fII.^2,ht_fIII.^2,ht_fIV.^2,ht_fV.^2,ht_fVI.^2,n,nn,'cor_int');

E=gp*totalh+hp*totalu;

end

