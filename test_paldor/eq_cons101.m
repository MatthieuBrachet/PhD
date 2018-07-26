function[eq2_fI, eq2_fII, eq2_fIII, eq2_fIV, eq2_fV, eq2_fVI]=eq_cons101(ht_fI, ht_fII,...
    ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
% conservation equation (right hand side) for the shallow water system.
global n nn
global hp

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div101(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

eq2_fI=-hp*div_fI;
eq2_fII=-hp*div_fII;
eq2_fIII=-hp*div_fIII;
eq2_fIV=-hp*div_fIV;
eq2_fV=-hp*div_fV;
eq2_fVI=-hp*div_fVI;