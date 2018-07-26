function [ M ] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI )
%mass
global nrm
global n nn

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,M]=nrm101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn,nrm);

end

