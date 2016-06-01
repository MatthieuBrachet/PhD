function [ Qx,Qy,Qz ] = qte_mvt(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn)
%% quantité en x
qtx_fI=ht_fI.*vt_fI(:,:,1);
qtx_fII=ht_fII.*vt_fII(:,:,1);
qtx_fIII=ht_fIII.*vt_fIII(:,:,1);
qtx_fIV=ht_fIV.*vt_fIV(:,:,1);
qtx_fV=ht_fV.*vt_fV(:,:,1);
qtx_fVI=ht_fVI.*vt_fVI(:,:,1);

[~,~,~,~,~,~,Qx]=...
    nrm72(qtx_fI,qtx_fII,qtx_fIII,qtx_fIV,qtx_fV,qtx_fVI,n,nn,'int');

%% quantité en y
qty_fI=ht_fI.*vt_fI(:,:,2);
qty_fII=ht_fII.*vt_fII(:,:,2);
qty_fIII=ht_fIII.*vt_fIII(:,:,2);
qty_fIV=ht_fIV.*vt_fIV(:,:,2);
qty_fV=ht_fV.*vt_fV(:,:,2);
qty_fVI=ht_fVI.*vt_fVI(:,:,2);

[~,~,~,~,~,~,Qy]=...
    nrm72(qty_fI,qty_fII,qty_fIII,qty_fIV,qty_fV,qty_fVI,n,nn,'int');

%% quantité en z
qtz_fI=ht_fI.*vt_fI(:,:,3);
qtz_fII=ht_fII.*vt_fII(:,:,3);
qtz_fIII=ht_fIII.*vt_fIII(:,:,3);
qtz_fIV=ht_fIV.*vt_fIV(:,:,3);
qtz_fV=ht_fV.*vt_fV(:,:,3);
qtz_fVI=ht_fVI.*vt_fVI(:,:,3);

[~,~,~,~,~,~,Qz]=...
    nrm72(qtz_fI,qtz_fII,qtz_fIII,qtz_fIV,qtz_fV,qtz_fVI,n,nn,'int');
end

