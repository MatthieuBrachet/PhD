close all
%% courbes

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72( v_fI  , v_fII  , v_fIII  , v_fIV  , v_fV  , v_fVI  , n , nn);
    
figure(8)
plot_cs11(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI);
%%
[divb_fI,divb_fII,divb_fIII,divb_fIV,divb_fV,divb_fVI]=...
        div72( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI  , n , nn);
    
figure(9)
plot_cs11(n,nn,divb_fI,divb_fII,divb_fIII,divb_fIV,divb_fV,divb_fVI);


%%
figure(11)
subplot(131)
plot_cs11(n,nn,vt_fI(:,:,1),vt_fII(:,:,1),vt_fIII(:,:,1),vt_fIV(:,:,1),vt_fV(:,:,1),vt_fVI(:,:,1));
title('x-component')

subplot(132)
plot_cs11(n,nn,vt_fI(:,:,2),vt_fII(:,:,2),vt_fIII(:,:,2),vt_fIV(:,:,2),vt_fV(:,:,2),vt_fVI(:,:,2));
title('y-component')

subplot(133)
plot_cs11(n,nn,vt_fI(:,:,3),vt_fII(:,:,3),vt_fIII(:,:,3),vt_fIV(:,:,3),vt_fV(:,:,3),vt_fVI(:,:,3));
title('z-component')