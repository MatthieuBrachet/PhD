function[eq2_fI, eq2_fII, eq2_fIII, eq2_fIV, eq2_fV, eq2_fVI]=eq_cons74(ht_fI, ht_fII,...
    ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)

global n nn

% face I
vec_I=zeros(size(vt_fI));
for i=1:3
    vec_I(:,:,i)=ht_fI(:,:).*vt_fI(:,:,i);
end

% face II
vec_II=zeros(size(vt_fII));
for i=1:3
    vec_II(:,:,i)=ht_fII(:,:).*vt_fII(:,:,i);
end
% face III
vec_III=zeros(size(vt_fIII));
for i=1:3
    vec_III(:,:,i)=ht_fIII(:,:).*vt_fIII(:,:,i);
end
% face IV
vec_IV=zeros(size(vt_fIV));
for i=1:3
    vec_IV(:,:,i)=ht_fIV(:,:).*vt_fIV(:,:,i);
end
% face V
vec_V=zeros(size(vt_fV));
for i=1:3
    vec_V(:,:,i)=ht_fV(:,:).*vt_fV(:,:,i);
end
% face VI
vec_VI=zeros(size(vt_fVI));
for i=1:3
    vec_VI(:,:,i)=ht_fVI(:,:).*vt_fVI(:,:,i);
end

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div74(vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI,n,nn);

eq2_fI=-div_fI;
eq2_fII=-div_fII;
eq2_fIII=-div_fIII;
eq2_fIV=-div_fIV;
eq2_fV=-div_fV;
eq2_fVI=-div_fVI;