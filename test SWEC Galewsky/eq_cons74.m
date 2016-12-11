function[eq2_fI, eq2_fII, eq2_fIII, eq2_fIV, eq2_fV, eq2_fVI]=eq_cons74(ht_fI, ht_fII,...
    ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)

global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
global visc

% face I
vec_I=zeros(size(vt_fI));
[hs_fI] = relief(x_fI,y_fI,z_fI);
hstar_fI=ht_fI-hs_fI;
for i=1:3
    vec_I(:,:,i)=hstar_fI(:,:).*vt_fI(:,:,i);
end

% face II
vec_II=zeros(size(vt_fII));
[hs_fII] = relief(x_fII,y_fII,z_fII);
hstar_fII=ht_fII-hs_fII;
for i=1:3
    vec_II(:,:,i)=hstar_fII(:,:).*vt_fII(:,:,i);
end
% face III
vec_III=zeros(size(vt_fIII));
[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
hstar_fIII=ht_fIII-hs_fIII;
for i=1:3
    vec_III(:,:,i)=hstar_fIII(:,:).*vt_fIII(:,:,i);
end
% face IV
vec_IV=zeros(size(vt_fIV));
[hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
hstar_fIV=ht_fIV-hs_fIV;
for i=1:3
    vec_IV(:,:,i)=hstar_fIV(:,:).*vt_fIV(:,:,i);
end
% face V
vec_V=zeros(size(vt_fV));
[hs_fV] = relief(x_fV,y_fV,z_fV);
hstar_fV=ht_fV-hs_fV;
for i=1:3
    vec_V(:,:,i)=hstar_fV(:,:).*vt_fV(:,:,i);
end
% face VI
vec_VI=zeros(size(vt_fVI));
[hs_fVI] = relief(x_fVI,y_fVI,z_fVI);
hstar_fVI=ht_fVI-hs_fVI;
for i=1:3
    vec_VI(:,:,i)=hstar_fVI(:,:).*vt_fVI(:,:,i);
end

[grl_fI,grl_fII,grl_fIII,grl_fIV,grl_fV,grl_fVI]=...
    gr74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn);
[lap_fI,lap_fII,lap_fIII,lap_fIV,lap_fV,lap_fVI]=...
    div74(grl_fI,grl_fII,grl_fIII,grl_fIV,grl_fV,grl_fVI,n,nn);

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div74(vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI,n,nn);

eq2_fI=-div_fI+visc*lap_fI;
eq2_fII=-div_fII+visc*lap_fII;
eq2_fIII=-div_fIII+visc*lap_fIII;
eq2_fIV=-div_fIV+visc*lap_fIV;
eq2_fV=-div_fV+visc*lap_fV;
eq2_fVI=-div_fVI+visc*lap_fVI;