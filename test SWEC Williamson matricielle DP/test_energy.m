% test energie locale

global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
global gp


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

[n1,n2]=size(x_fI);
for i=1:n1
    for j=1:n2
        norm_I(i,j)=dot(hstar_fI(i,j)*vt_fI(i,j,:),vt_fI(i,j,:));
        norm_II(i,j)=dot(hstar_fII(i,j)*vt_fII(i,j,:),vt_fII(i,j,:));
        norm_III(i,j)=dot(hstar_fIII(i,j)*vt_fIII(i,j,:),vt_fIII(i,j,:));
        norm_IV(i,j)=dot(hstar_fIV(i,j)*vt_fIV(i,j,:),vt_fIV(i,j,:));
        norm_V(i,j)=dot(hstar_fV(i,j)*vt_fV(i,j,:),vt_fV(i,j,:));
        norm_VI(i,j)=dot(hstar_fVI(i,j)*vt_fVI(i,j,:),vt_fVI(i,j,:));
    end
end


en_I=.5*norm_I+.5*gp*(ht_fI.^2-hs_fI.^2);
en_II=.5*norm_II+.5*gp*(ht_fII.^2-hs_fII.^2);
en_III=.5*norm_III+.5*gp*(ht_fIII.^2-hs_fIII.^2);
en_IV=.5*norm_IV+.5*gp*(ht_fIV.^2-hs_fIV.^2);
en_V=.5*norm_V+.5*gp*(ht_fV.^2-hs_fV.^2);
en_VI=.5*norm_VI+.5*gp*(ht_fVI.^2-hs_fVI.^2);


hFig = figure(20);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs102(n,nn,log10(abs(en_I)+eps),log10(abs(en_II)+eps),log10(abs(en_III)+eps),log10(abs(en_IV)+eps),log10(abs(en_V)+eps),log10(abs(en_VI)+eps))
title(['log10(abs(energy)) at time : ', num2str(time(end))])
colorbar


hFig = figure(21);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs102(n,nn,0.5*norm_I,.5*norm_II,.5*norm_III,.5*norm_IV,.5*norm_V,.5*norm_VI)
title(['energy kinetic at time : ', num2str(time(end))])
colorbar

fig_placier
