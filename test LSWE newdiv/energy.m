function [ E ] = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn)
global gp hp

[pp,qq,~]=size(vt_fI);
for i=1:pp
    for j=1:qq
        norm_I(i,j)=dot(vt_fI(i,j,:),vt_fI(i,j,:));
        norm_II(i,j)=dot(vt_fII(i,j,:),vt_fII(i,j,:));
        norm_III(i,j)=dot(vt_fIII(i,j,:),vt_fIII(i,j,:));
        norm_IV(i,j)=dot(vt_fIV(i,j,:),vt_fIV(i,j,:));
        norm_V(i,j)=dot(vt_fV(i,j,:),vt_fV(i,j,:));
        norm_VI(i,j)=dot(vt_fVI(i,j,:),vt_fVI(i,j,:));
    end
end

[~,~,~,~,~,~,totalu]=...
    nrm72(norm_I,norm_II,norm_III,norm_IV,norm_V,norm_VI,n,nn,'int');

[~,~,~,~,~,~,totalh]=...
    nrm72(ht_fI.^2,ht_fII.^2,ht_fIII.^2,ht_fIV.^2,ht_fV.^2,ht_fVI.^2,n,nn,'int');

E=gp*totalh+hp*totalu;

end

